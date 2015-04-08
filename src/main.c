/*
 * HTseqQA
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
KSEQ_INIT(gzFile,gzread);


// pre declare function 
void generate_report();
void generate_report_gs();
int gc_content(std::string);


// Our arguments struct
struct args_t{
	char* incoming;
	int offset = -1;
	int max_quality_score = -1;
	/*outputs*/
	int report_only = 0;
	int greyscale = 0;
	int force = 0;
	/*file paths*/
	std::string report; //yaml report
	std::string R; // R report
	std::string basename; // basename + no suffix
	std::string no_suffix;

	}args;





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



// Strign for opt_arg
static const char *optString = "i:ho:rgf";

// Usage statement//
void usage(){
	std::cout << "./HTseqQA -i <fastq> " << std::endl;
	std::cout << "Options : " << std::endl;
	std::cout << "\t -o <offset> , manually set your offset" << std::endl;
	std::cout << "\t -r , Only print out the Rscript" << std::endl;
	std::cout << "\t -g , Set graphs for grey scale only" << std::endl;
	std::cout << "\t -f , Force the re run if tables already exist" << std::endl;
	return;
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
		case 'r':
			args.report_only = 1;
			break;
		case 'g':
			args.greyscale = 1;
			break;
		case 'f':
			args.force = 1;
			break;
		default:
			usage();
			exit(200);
		}

	opt = getopt(argc,argv,optString);
	}	

if (args.force == 0){
	/*try to open a file*/
	std::string file = args.no_suffix+"_gcdistributions.tbl";
	std::ifstream lastOutput(file);
	if (lastOutput){
		std::cerr<< args.basename << " appears to have been processed already, use -f to re-write tables" <<std::endl;
		exit(1);
		}
	}

if (args.report_only > 0){
	if (args.greyscale){
		generate_report_gs();
		}
	else{
		generate_report();
		}
	exit(1);
	}
// Provided arguments //1
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
std::unordered_map<int,unsigned long long> GC_distributions;

std::cout << args.basename << " started!" << std::endl;
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
 	GC_distributions[gc_content(sequence)]++;
 	/*Passing_Reads Processor*/
 	/*ascending sort*/
 	std::sort(quals.begin(),quals.end());

 	for (int i = 0; i < end; ++i){
 		for (int j = 0; j <= ((int)quals[i]-args.offset); ++j){
 			Passing_Reads[end-i-1][j]++;
 			}
 		}
	}

// Done with our Fastq file

/* Get our largest seen quality score*/
for(auto it : Quals_Distribution) {
   	 if (it.first > args.max_quality_score){
   	  args.max_quality_score = it.first;	
   	 }  
	} 

// How many basepairs did I see ?
long double total_seen = 0;
for (int i = 0; i < args.max_quality_score+1; ++i){ total_seen+=Quals_Distribution[i];}

// Print out how many reads we saw
std::ofstream reads(args.no_suffix+"_readcount.txt");
reads << read_counter << std::endl;
reads.close();

//print out my Cumulative qualities 
std::ofstream cum_quals(args.no_suffix+"_cum_quals.tbl");
cum_quals << "label\txlabel\tylabel"<<std::endl;
unsigned long long current_seen = 0;
for ( int j = 0; j <=args.max_quality_score; ++j){
	current_seen += Quals_Distribution[j];
	cum_quals << args.basename <<"\t"<< j << "\t"<< (current_seen/total_seen) << std::endl;
	} 
cum_quals.close();

//print out our nucleotides by position table
std::string Nucleotides ("ATGCN");
std::ofstream stacked_bargraph(args.no_suffix+"_nucleotides_by_position.tbl");
int read_length = NTs_ByPosition.size();
stacked_bargraph << "position\tnucleotide\tvalue\tlabel"<<std::endl;
for (int i = 0; i < read_length; ++i){
	for (int j =0; j<5;++j){
		stacked_bargraph << i+1 << "\t" <<Nucleotides[j]<<"\t"<< (NTs_ByPosition[i][Nucleotides[j]]/ read_counter ) << "\t"<< args.basename << std::endl;
		}
	}
stacked_bargraph.close();


//Print out our quality distributions
std::ofstream violin_plot(args.no_suffix+"_quality_distribution.tbl");
int quals_length = Quals_ByPosition.size();
violin_plot << "xlabel\tylabel\ty_weight\tlabel"<<std::endl;
for (int i = 0; i < quals_length ; ++i){	
	for (int j = 0; j < args.max_quality_score+1; ++j){
		violin_plot << i << "\t" << j <<"\t" << (Quals_ByPosition[i][j] / read_counter )<<"\t"<<args.basename<<std::endl;
		}
	}
violin_plot.close();

//Print out our passing reads filter
std::ofstream passing_reads(args.no_suffix+"_passing_reads.tbl");
int passing_reads_size = Passing_Reads.size();
passing_reads<<"xlabel\tylabel\tvalue" << std::endl;
for (int x =0; x < passing_reads_size; ++x){
	for (int y = 0; y <= args.max_quality_score; ++y){
		if ((Passing_Reads[x][y]/read_counter) == 0){
			passing_reads<< x+1 << "\t" << y <<"\tNA"<<std::endl;
			}
		else{
			passing_reads<< x+1 << "\t" << y <<"\t"<<(Passing_Reads[x][y]/read_counter)<<std::endl;
			}
	
		}
	}
passing_reads.close();

// Print out our novelty table
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

//print out our GC distributions table
std::ofstream gc_distribution(args.no_suffix+"_gcdistributions.tbl");
gc_distribution<<"GC\tcount"<<std::endl;
for (int i = 0; i < 101 ; ++i){
	gc_distribution << i << "\t" << GC_distributions[i] << std::endl;
	}
gc_distribution.close();

if (args.greyscale){
	generate_report_gs();
	}
else{
	generate_report();
	}

return 0;
}


int gc_content(std::string incoming){
	int length = incoming.size();
	int gc = 0;
	for (int i = 0; i < length;++i ){
		if (incoming[i] == 'G' || incoming[i] == 'C'){
			gc++;
			}
		}
	return int(gc*100/length);
	}

void generate_report_gs(){

	std::ofstream R(args.R);

	std::string cum_quals_file_name = std::string(args.basename)+"_cum_quals.tbl"; 
	std::string stacked_bargraph_file_name = std::string(args.basename)+"_nucleotides_by_position.tbl";
	std::string violin_file_name = std::string(args.basename)+"_quality_distribution.tbl";
	std::string passing_reads_file_name = std::string(args.basename)+"_passing_reads.tbl"; 
	std::string novelty_file_name = std::string(args.basename)+"_novelty.tbl";
	std::string gcdistribution_file_name = std::string(args.basename)+"_gcdistributions.tbl";

	R << "library(ggplot2)" << std::endl;

	R << "df <- as.data.frame(read.table(\""<<cum_quals_file_name<<"\",header=TRUE))"<<std::endl;
	R << "p <-ggplot(df, aes (x = xlabel, y=ylabel,color =label,group =label))";
	R << "+ geom_point(color=\"black\")+geom_line(color=\"black\") + ";
	R << "theme(axis.text.x = element_text(size = 3 ,angle=90),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
	R <<"+ ggtitle(\"Cumulative proportion of qualities" << std::endl;
	R << "for " << std::string(args.basename) <<"\")"<< +"+xlab(\"Quality Score\")+ylab(\"Cumulative Proportion\") + labs(color = \"sample\")+ theme(axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))"<<std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".cqs.ps\",dpi=600,width=7,height=5)"<<std::endl;



	R << "df2 <- as.data.frame(read.table(\""<<stacked_bargraph_file_name<<"\",header=TRUE))"<<std::endl;
	R << "p <- ggplot(df2, aes (x = position, y=value ,fill = nucleotide),position=dodge)";
	R << "+geom_bar(stat=\"identity\",linetype=\"blank\")+scale_fill_grey(start=0,end=0.9)";
	R << "+theme(axis.text.x = element_text(angle=90),plot.background = element_blank(),panel.grid.major = element_blank()";
	R << ",panel.grid.minor = element_blank(),panel.background = element_blank()) ";
	R << "+ ggtitle(\"Nucleotide Distribution for "<< args.basename <<"\") ";
	R << " +xlab(\"Position\")+ylab(\"Proportion\")+labs(fill = \"Nucleotide Class\")+ theme(axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))"	 << std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".ntbps.ps\",dpi=600,width=7,height=5)"<<std::endl;



	R << "df3 <- as.data.frame(read.table(\""<<violin_file_name<<"\",header=TRUE))"<<std::endl;
	R <<  "p <- ggplot(df3, aes(x= factor(xlabel),y=ylabel,fill=label, weight = y_weight)) ";
	R << "+ geom_violin(scale=\"area\",position=\"identity\",linetype=\"blank\",fill=\"black\",colour=\"white\") ";
	R << "+theme(axis.text.x = element_text(size = 3, angle=90),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
	R << "+ ggtitle(\"Quality Score Distributions by Position\")+ xlab(\"Position\") ";
	R << "+ ylab(\"Quality Scores\") + labs(\""<< args.basename <<"\")+ theme(axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))" << std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".qdbs.ps\",dpi=600,width=7,height=5)"<<std::endl;


	R << "df4 <- as.data.frame(read.table(\""<<passing_reads_file_name<<"\",header=TRUE))"<<std::endl;
	R << "p <- ggplot(df4, aes(x = xlabel, y = ylabel , fill = value)) + geom_tile(aes(fill = value)) ";
	R << " + scale_fill_gradient(low=\"dark grey\", high=\"white\", na.value=\"black\")+theme(axis.text.x = element_text(angle=90),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
	R << " + ggtitle(\"Passing reads for "<<args.basename<<"\") ";
	R << "  + ylab(\"Minimum Quality Score\")+xlab(\"Number of bases >= Minimum Quality\")";
	R << " +labs(fill = \"Proportion of Reads\")+ theme(axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))" << std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".prf.ps\",dpi=600,width=7,height=5)"<<std::endl;

	R << "df5 <- as.data.frame(read.table(\""<<novelty_file_name<<"\",header=TRUE))"<<std::endl;
	R << "p <- ggplot(df5, aes(x=novelty,y=count),position=dodge)+geom_histogram(stat=\"identity\")+ggtitle(\"Uniqueness of reads for "<< args.basename << "\")";
	R << "+ylab (\"Number of reads in category\")+xlab(\"Uniqueness score\")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.background= element_blank(),axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))"<<std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".novelty.ps\",dpi=600,width=7,height=5)"<<std::endl;

	R << "df6 <- as.data.frame(read.table(\""<<gcdistribution_file_name<<"\",header=TRUE))"<<std::endl;
	R << "p <- ggplot(df6, aes(x=GC,y=count),position=dodge)+geom_histogram(stat=\"identity\")+ggtitle(\"Read GC content distribution for "<< args.basename << "\")";
	R << "+ylab (\"Number of reads\")+xlab(\"Read GC content\")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.background= element_blank(),axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))"<<std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".gcdist.ps\",dpi=600,width=7,height=5)"<<std::endl;


	R.close();

	std::cout << args.basename << " completed!" << std::endl;

	}

void generate_report(){

	std::ofstream R(args.R);

	std::string cum_quals_file_name = std::string(args.basename)+"_cum_quals.tbl"; 
	std::string stacked_bargraph_file_name = std::string(args.basename)+"_nucleotides_by_position.tbl";
	std::string violin_file_name = std::string(args.basename)+"_quality_distribution.tbl";
	std::string passing_reads_file_name = std::string(args.basename)+"_passing_reads.tbl"; 
	std::string novelty_file_name = std::string(args.basename)+"_novelty.tbl";
	std::string gcdistribution_file_name = std::string(args.basename)+"_gcdistributions.tbl";

	R << "library(ggplot2)" << std::endl;

	R << "df <- as.data.frame(read.table(\""<<cum_quals_file_name<<"\",header=TRUE))"<<std::endl;
	R << "p <-ggplot(df, aes (x = xlabel, y=ylabel,color =label,group =label))";
	R << "+ geom_point()+geom_line() + ";
	R << "theme(axis.text.x = element_text(size = 3 ,angle=90),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
	R <<"+ ggtitle(\"Cumulative proportion of qualities" << std::endl;
	R << "for " << std::string(args.basename) <<"\")"<< +"+xlab(\"Quality Score\")+ylab(\"Cumulative Proportion\") + labs(color = \"sample\")+ theme(axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))"<<std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".cqs.ps\",dpi=600,width=7,height=5)"<<std::endl;



	R << "df2 <- as.data.frame(read.table(\""<<stacked_bargraph_file_name<<"\",header=TRUE))"<<std::endl;
	R << "p <- ggplot(df2, aes (x = position, y=value ,fill = nucleotide),position=dodge)";
	R << "+geom_bar(stat=\"identity\",linetype=\"blank\")";
	R << "+theme(axis.text.x = element_text(angle=90),plot.background = element_blank(),panel.grid.major = element_blank()";
	R << ",panel.grid.minor = element_blank(),panel.background = element_blank()) ";
	R << "+ ggtitle(\"Nucleotide Distribution for "<< args.basename <<"\") ";
	R << " +xlab(\"Position\")+ylab(\"Proportion\")+labs(fill = \"Nucleotide Class\")+ theme(axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))"	 << std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".ntbps.ps\",dpi=600,width=7,height=5)"<<std::endl;



	R << "df3 <- as.data.frame(read.table(\""<<violin_file_name<<"\",header=TRUE))"<<std::endl;
	R <<  "p <- ggplot(df3, aes(x= factor(xlabel),y=ylabel,fill=label, weight = y_weight)) ";
	R << "+ geom_violin(scale=\"area\",position=\"identity\",linetype=\"blank\") ";
	R << "+theme(axis.text.x = element_text(size = 3, angle=90),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
	R << "+ ggtitle(\"Quality Score Distributions by Position\")+ xlab(\"Position\") ";
	R << "+ ylab(\"Quality Scores\") + labs(\""<< args.basename <<"\")+ theme(axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))" << std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".qdbs.ps\",dpi=600,width=7,height=5)"<<std::endl;


	R << "df4 <- as.data.frame(read.table(\""<<passing_reads_file_name<<"\",header=TRUE))"<<std::endl;
	R << "p <- ggplot(df4, aes(x = xlabel, y = ylabel , fill = value)) + geom_tile(aes(fill = value)) ";
	R << " + scale_fill_gradient(low=\"blue\", high=\"white\")+theme(axis.text.x = element_text(angle=90),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
	R << " + ggtitle(\"Passing reads for "<<args.basename<<"\") ";
	R << "  + ylab(\"Minimum Quality Score\")+xlab(\"Number of bases >= Minimum Quality\")";
	R << " +labs(fill = \"Proportion of Reads\")+ theme(axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))" << std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".prf.ps\",dpi=600,width=7,height=5)"<<std::endl;

	R << "df5 <- as.data.frame(read.table(\""<<novelty_file_name<<"\",header=TRUE))"<<std::endl;
	R << "p <- ggplot(df5, aes(x=novelty,y=count),position=dodge)+geom_histogram(stat=\"identity\")+ggtitle(\"Uniqueness of reads for "<< args.basename << "\")";
	R << "+ylab (\"Number of reads in category\")+xlab(\"Uniqueness score\")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.background= element_blank(),axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))"<<std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".novelty.ps\",dpi=600,width=7,height=5)"<<std::endl;

	R << "df6 <- as.data.frame(read.table(\""<<gcdistribution_file_name<<"\",header=TRUE))"<<std::endl;
	R << "p <- ggplot(df6, aes(x=GC,y=count),position=dodge)+geom_histogram(stat=\"identity\")+ggtitle(\"Read GC content distribution for "<< args.basename << "\")";
	R << "+ylab (\"Number of reads\")+xlab(\"Read GC content\")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.background= element_blank(),axis.text.x = element_text(colour=\"black\"), axis.text.y = element_text(colour=\"black\") , axis.title.x= element_text(vjust=-0.5) , axis.title.y= element_text(vjust=0.25))"<<std::endl;
	R << "ggsave(p,file=\"" << args.basename <<".gcdist.ps\",dpi=600,width=7,height=5)"<<std::endl;


	R.close();

	std::cout << args.basename << " completed!" << std::endl;

	}