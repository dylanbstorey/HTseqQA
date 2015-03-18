/*
 * untitled.c
 * 
 * Copyright 2014 Dylan Storey <dstorey@optimus>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <stdio.h>
#include <kseq.h>
#include <zlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <set>
#include <libgen.h>
#include <algorithm>


// Initialize for kseq
KSEQ_INIT(gzFile,gzread)

// Our arguments struct
struct args_t{
	char* incoming;
	int offset = -1;
	int max_quality_score = -1;
	/*outputs*/

	/*file paths*/
	std::string report; //yaml report
	std::string R; // R report
	std::string basename; // basename + no suffix
	std::string no_suffix;
	

	}args;

// Strign for opt_arg
static const char *optString = "i:ho:";

// Usage statement//
void usage(){
	std::cout << "./fastqcreport -i fastq" << std::endl;
	return;
	}

// Function for guessing a files offset
int get_offset(char* file){
	/*open file*/
	gzFile inc = gzopen(file , "r");
	kseq_t *records = kseq_init(inc);
	int l;
	while ((l = kseq_read(records)) >= 0 ){
		std::string quals = records->qual.s;
		for(std::string::size_type i = 0; i < quals.size(); ++i){
			if (quals[i] < 58){
				kseq_destroy(records);
				gzclose(inc);
				return 33;
				}
			else if (quals[i] > 74){
				kseq_destroy(records);
				gzclose(inc);
				return 64;
				}
			}
		}
	kseq_destroy(records);
	gzclose(inc);
	return -1;
	}




int main(int argc, char **argv){

// Get options
int opt = getopt(argc,argv,optString);

while (opt != -1){
	switch (opt){
		case 'i':
			args.incoming = optarg;
			
			args.no_suffix = optarg;
			args.no_suffix = args.no_suffix.substr(0,args.no_suffix.find_first_of("."));
	
			args.R = args.no_suffix + ".R"; // R report
			
			args.basename = args.incoming; // basename + no suffix
			args.basename = basename(args.incoming);
			args.basename = args.basename.substr(0,args.basename.find_first_of("."));

			break;
		case 'h':
			usage();
			exit(0);
		case 'o':
			args.offset = atoi(optarg);
			break;
		default:
			usage();
			exit(200);
		}

	opt = getopt(argc,argv,optString);
	}	


// Retard Checking //1
if (not args.incoming){
	std::cerr << "No fastq provided !" << std::endl;
	exit(200);
	}
// Quality Offset Checking
if (args.offset == -1){
	args.offset = get_offset(args.incoming);
	if (args.offset == -1){
		std::cerr << "Couldn't guess your offset! use the -o option to set manually" << std::endl;
		exit(200);
		}
	}



gzFile inc = gzopen(args.incoming , "r");
kseq_t *records = kseq_init(inc);
int l;


/* Containers for tracking data */
std::unordered_map < int, std::unordered_map<char , unsigned long long>> NTs_ByPosition;
std::unordered_map <int, std::unordered_map<int , unsigned long long>> Quals_ByPosition;
std::unordered_map<int, unsigned long long> Quals_Distribution;
long double read_counter = 0;
std::unordered_map<int, std::unordered_map<int,unsigned long long>> Passing_Reads;
std::unordered_map<std::string,int> Sequence_Representation;

while ((l = kseq_read(records)) >=0){
 	// we have records;

 	std::string header = records->name.s;
 	std::string sequence = records-> seq.s;
 	std::string quals = records->qual.s;
 	Sequence_Representation[sequence]++;
 	/*process reads*/
 	read_counter++;
 	int end = sequence.size();
 	for(int i = 0; i < end; ++i){
 		NTs_ByPosition[i][sequence[i]]++;
 		Quals_ByPosition[i][quals[i]-args.offset]++;
 		Quals_Distribution[quals[i]-args.offset]++;
 		}

 	/*Passing_Reads Processor*/
 	/*ascending sort*/
 	std::sort(quals.begin(),quals.end());

 	for (int i = 0; i < end; ++i){
 		for (int j = 0; j <= ((int)quals[i]-args.offset); ++j){
 			Passing_Reads[end-i-1][j]++;
 			}
 		}
	}

/*All reads are in*/

/* Get our largest seen quality score*/
for(auto it : Quals_Distribution) {
   	 if (it.first > args.max_quality_score){
   	  args.max_quality_score = it.first;	
   	 }  
	} 

// How many basepairs did I see ?
long double total_seen = 0;
for (int i = 0; i < args.max_quality_score+1; ++i){ total_seen+=Quals_Distribution[i];}

// My Nucleotide classes
std::string Nucleotides ("ATGCN");
std::ofstream cum_quals(args.no_suffix+"_cum_quals.tbl");
cum_quals << "label\txlabel\tylabel"<<std::endl;

unsigned long long current_seen = 0;
for ( int j = 0; j <=args.max_quality_score; ++j){
	current_seen += Quals_Distribution[j];
	cum_quals << args.basename <<"\t"<< j << "\t"<< (current_seen/total_seen) << std::endl;
	} 
cum_quals.close();

std::ofstream stacked_bargraph(args.no_suffix+"_nucleotides_by_position.tbl");
int read_length = NTs_ByPosition.size();
stacked_bargraph << "position\tnucleotide\tvalue\tlabel"<<std::endl;
for (int i = 0; i < read_length; ++i){
	for (int j =0; j<5;++j){
		stacked_bargraph << i+1 << "\t" <<Nucleotides[j]<<"\t"<< (NTs_ByPosition[i][Nucleotides[j]]/ read_counter ) << "\t"<< args.basename << std::endl;
		}
	}
stacked_bargraph.close();

std::ofstream violin_plot(args.no_suffix+"_quality_distribution.tbl");
int quals_length = Quals_ByPosition.size();
violin_plot << "xlabel\tylabel\ty_weight\tlabel"<<std::endl;
for (int i = 0; i < quals_length ; ++i){	
	
	for (int j = 0; j < args.max_quality_score+1; ++j){
		violin_plot << i << "\t" << j <<"\t" << (Quals_ByPosition[i][j] / read_counter )<<"\t"<<args.basename<<std::endl;
		}
	}

violin_plot.close();

std::ofstream passing_reads(args.no_suffix+"_passing_reads.tbl");
int passing_reads_size = Passing_Reads.size();

passing_reads<<"xlabel\tylabel\tvalue" << std::endl;


for (int x =0; x < passing_reads_size; ++x){
	for (int y = 0; y <= args.max_quality_score; ++y){
		if ((Passing_Reads[x][y]/read_counter) == 0){
			passing_reads<< x+1 << "\t" << y <<"\t-1"<<std::endl;
			}
		else{
			passing_reads<< x+1 << "\t" << y <<"\t"<<(Passing_Reads[x][y]/read_counter)<<std::endl;
			}
	
		}
	}
passing_reads.close();

std::ofstream novelty(args.no_suffix+"_novelty.tbl");
std::unordered_map <int , int> novelty_box;
for (auto sequence : Sequence_Representation){
	novelty_box[sequence.second]++;
	}

novelty << "novelty\tcount" << std::endl;

for (auto bin :  novelty_box){
	novelty << bin.first << "\t"<< bin.second << std::endl;
	}

novelty.close();




std::ofstream R(args.R);

std::string cum_quals_file_name = std::string(args.basename)+"_cum_quals.tbl"; 
std::string stacked_bargraph_file_name = std::string(args.basename)+"_nucleotides_by_position.tbl";
std::string violin_file_name = std::string(args.basename)+"_quality_distribution.tbl";
std::string passing_reads_file_name = std::string(args.basename)+"_passing_reads.tbl"; 
std::string novelty_file_name = std::string(args.basename)+"_novelty.tbl";

R << "library(ggplot2)" << std::endl;
R << "postscript(\""<< args.basename <<".cqs.ps\", width=2000, height=2000)" << std::endl;
R << "df <- as.data.frame(read.table(\""<<cum_quals_file_name<<"\",header=TRUE))"<<std::endl;
R << "ggplot(df, aes (x = xlabel, y=ylabel,color =label,group =label))";
R << "+ geom_point()+geom_line() + ";
R << "theme(axis.text.x = element_text(size=4,angle=90),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
R <<"+ ggtitle(\"Cumulative proportion of qualities" << std::endl;
R << "for " << std::string(args.basename) <<"\")"<< +"+xlab(\"Quality Score\")+ylab(\"Cumulative Proportion\") + labs(color = \"sample\")"<<std::endl;



R << "postscript(\""<< args.basename <<".ntbps.ps\", width=2000, height=2000)" << std::endl;
R << "df2 <- as.data.frame(read.table(\""<<stacked_bargraph_file_name<<"\",header=TRUE))"<<std::endl;
R << "ggplot(df2, aes (x = position, y=value ,fill = nucleotide))";
R << "+geom_bar(stat=\"identity\",linetype=\"blank\")";
R << "+theme(axis.text.x = element_text(size=4,angle=90),plot.background = element_blank(),panel.grid.major = element_blank()";
R << ",panel.grid.minor = element_blank(),panel.background = element_blank()) ";
R << "+ ggtitle(\"Nucleotide Distribution for "<< args.basename <<"\") ";
R << " +xlab(\"Position\")+ylab(\"Proportion\")+labs(fill = \"Nucleotide Class\")"	 << std::endl;

R << "postscript(\""<< args.basename <<".qdbs.ps\", width=2000, height=2000)" << std::endl;
R << "df3 <- as.data.frame(read.table(\""<<violin_file_name<<"\",header=TRUE))"<<std::endl;
R <<  "ggplot(df3, aes(x= factor(xlabel),y=ylabel,fill=label, weight = y_weight)) ";
R << "+ geom_violin(scale=\"count\",position=\"identity\",linetype=\"blank\") ";
R << "+theme(axis.text.x = element_text(size=4,angle=90),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
R << "+ ggtitle(\"Quality Score Distributions by Position\")+ xlab(\"Position\") ";
R << "+ ylab(\"Quality Scores\") + labs(\""<< args.basename <<"\")" << std::endl;

R << "postscript(\""<< args.basename <<".prf.ps\", width=2000, height=2000)" << std::endl;
R << "df4 <- as.data.frame(read.table(\""<<passing_reads_file_name<<"\",header=TRUE))"<<std::endl;
R << "ggplot(df4, aes(x = xlabel, y = ylabel , fill = value)) + geom_tile(aes(fill = value)) ";
R << " + scale_fill_gradient(low=\"steelblue\", high=\"white\")+theme(axis.text.x = element_text(size=4,angle=90),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
R << " + ggtitle(\"Passing reads for "<<args.basename<<"\") ";
R << "  + ylab(\"Minimum Quality Score\")+xlab(\"Number of bases >= Minimum Quality\")";
R << " +labs(fill = \"Proportion of Reads\")" << std::endl;

R << "postscript(\""<< args.basename <<".novelty.ps\", width=2000, height=2000)" << std::endl;
R << "df5 <- as.data.frame(read.table(\""<<novelty_file_name<<"\",header=TRUE))"<<std::endl;
R << "ggplot(df5, aes(x=novelty,y=count))+geom_histogram(stat=\"identity\")+ggtitle(\"Uniqueness of reads for "<< args.basename << "\")";
R << "+ylab (\"Number of reads in category\")+xlab(\"Uniqueness score\")"<<std::endl;


R.close();








// R << "df4 <- data.frame (" << x_axis << " , "<< y_axis << " , "<< value <<")"<< std::endl;


return 0;
}
