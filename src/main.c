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
	int offset = 0;
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
			if (quals[i] < 59){
				kseq_destroy(records);
				gzclose(inc);
				return 33;
				}
			else if (quals[i] > 73){
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
			
			args.report = args.no_suffix + ".report.yaml";

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
if (not args.offset){
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


while ((l = kseq_read(records)) >=0){
 	// we have records;

 	std::string header = records->name.s;
 	std::string sequence = records-> seq.s;
 	std::string quals = records->qual.s;
 	
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
 			Passing_Reads[end-i][j]++;
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

// /*Dump to a YAML*/
// std::ofstream ReportYAML(args.report);
// ReportYAML << "---" << std::endl;
// ReportYAML << "NTs_by_position:"<<std::endl;

// /* Dump our NTs_By_Position structure */
// std::string Nucleotides ("ATGCN");
// for (int i = 0; i < NTs_ByPosition.size() ; ++i){
		
// 		for (int j =0; j < Nucleotides.size(); ++j){
// 			if (j == 0){
// 				ReportYAML << "- - ";
// 				}
// 			else{
// 				ReportYAML << "  - ";
// 				}
// 			ReportYAML << (double(NTs_ByPosition[i][Nucleotides[j]])/ read_counter ) << std::endl;
//  		}

// 	}

// /*Dump total number of reads*/
// ReportYAML << "number_of_reads: " << read_counter << std::endl;

// /*Dump cumulative quality scores*/
// ReportYAML << "cumulative_quals:" << std::endl;
// unsigned long current_seen = 0;
// float total_seen = 0;
// for (int i = 0; i < args.max_quality_score+1; ++i){ total_seen+=Quals_Distribution[i];}
// for (int i = 0; i < args.max_quality_score+1 ; ++i){
// 	current_seen += Quals_Distribution[i];
// 	ReportYAML << "- " << double(current_seen)/total_seen << std::endl;
// 	}


// /*Dump our quality scores by position*/
// ReportYAML<<"quals_by_position:" << std::endl;
// for (int i = 0; i < Quals_ByPosition.size() ; ++i){	
// 		for (int j = 0; j < args.max_quality_score+1; ++j){
// 			if (j == 0){
// 				ReportYAML << "- - ";
// 				}
// 			else{
// 				ReportYAML << "  - ";
// 				}
// 			ReportYAML << (	double(Quals_ByPosition[i][j]) / read_counter ) << std::endl;
//  		}
// 	}

// /*Dump our passing reads filter*/

// ReportYAML<<"passing_reads:"<<std::endl;
// for (int i =0; i < Passing_Reads.size(); ++i){
// 	ReportYAML << "- ";
// 	for (int j = 0; j <= args.max_quality_score; ++j){
// 		ReportYAML << "- " << (double(Passing_Reads[i][j])/ read_counter) << std::endl ;
// 		}
// 	}

// ReportYAML.close();

// */

/*
Rfile script outputs
*/

/*cumulative qualitites*/
std::ofstream R(args.R);
R << "library(ggplot2)" << std::endl;
R << "postscript(\""<< args.basename <<".qc.ps\", width=500, height=500)" << std::endl;
R << "df <- data.frame(labels = c(" ;
/*print out our labels*/
std::string name = args.basename;
name ="\"" + name + "\"";

for (int j =0; j <= args.max_quality_score; ++j){
	if (j < args.max_quality_score){
		R<< name << ",";
		}
	else {
		R << name << "),";
		}
	}
	
/*print out our x_values*/
R << "x_label = factor(c(";
for (int j = 0; j <= args.max_quality_score ; ++j){
	if (j < args.max_quality_score){
		R << j << ",";
		}
	else{
		R << j<<")),";
		}
	}	
/*print out our y_values*/
R << "y_label = c(";
unsigned long long current_seen = 0;
for (int j = 0; j <= args.max_quality_score ; ++j){
	current_seen += Quals_Distribution[j];
	if (j < args.max_quality_score){
		R << (current_seen/total_seen) << ",";
		}
	else{
		R << (current_seen/total_seen) << "))" << std::endl;
		}
	}
R << "ggplot(df, aes (x = x_label, y=y_label,color =labels,group =labels))";
R << "+ geom_point()+geom_line() + ";
R << "theme(axis.text.x = element_text(size=4,angle=90),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
R <<"+ ggtitle(\"Cumulative proportion of qualities" << std::endl;
R << "for " << std::string(args.basename) <<"\")"<< +"+xlab(\"Quality Score\")+ylab(\"Cumulative Proportion\") + labs(color = \"sample\")"<<std::endl;

/*Stacked Bargraph*/
std::string labels = "labels = c(";
std::string enum_lab = "enum_lab=c(";
std::string val_lab = "val_lab =c(";

int read_lenght = NTs_ByPosition.size();

for (int i = 0; i < read_lenght ; ++i){
		/*val_lab = val_lab , Nucleotides[j] = enum_lab name = labels*/		
		for (int j =0; j < 5; ++j){
			if (i == read_lenght-1){
				
				if (Nucleotides[j] == 'N'){
					enum_lab +="\"";
					enum_lab += Nucleotides[j];
					enum_lab +="\"),";

					labels += std::to_string(i+1);
					labels += "),";

					val_lab += std::to_string(NTs_ByPosition[i][Nucleotides[j]]/ read_counter ) +")";
					}
				else{
					enum_lab += "\"";
					enum_lab +=Nucleotides[j];
					enum_lab += +"\",";

					labels += std::to_string(i+1); 
					labels += ",";

					val_lab += std::to_string(NTs_ByPosition[i][Nucleotides[j]]/ read_counter ) +",";
					}

				}
			else{
				labels += std::to_string(i+1);
				labels += ",";

				enum_lab += "\"";
				enum_lab += Nucleotides[j];
				enum_lab +="\",";
				
				val_lab += std::to_string(NTs_ByPosition[i][Nucleotides[j]]/ read_counter ) +",";
				}
 		}

	}

R<<"df2 <- data.frame(" << labels << enum_lab << val_lab <<")"<<std::endl;
R << "ggplot(df2, aes (x = labels, y=val_lab ,fill = enum_lab))";
R << "+geom_bar(stat=\"identity\",linetype=\"blank\")";
R << "+theme(axis.text.x = element_text(size=4,angle=90),plot.background = element_blank(),panel.grid.major = element_blank()";
R << ",panel.grid.minor = element_blank(),panel.background = element_blank()) ";
R << "+ ggtitle(\"Nucleotide Distribution for "<< args.basename <<"\") ";
R << " +xlab(\"Position\")+ylab(\"Proportion\")+labs(fill = \"Nucleotide Class\")"	 << std::endl;



/*violin plot qualities by position*/
std::string y_axis  = "y_label = c(";
std::string x_axis  = "x_label = factor(c(";
std::string y_value = "y_weights = c(";
std::string label   = "labels = c(";

int quals_length = Quals_ByPosition.size();
/*data fram expansion*/
for (int i = 0; i < quals_length ; ++i){	
		for (int j = 0; j < args.max_quality_score+1; ++j){
			if (i  == quals_length -1 ){
				if (j == args.max_quality_score){
					label += name;
					label += ")";
					
					x_axis += std::to_string(i+1);
					x_axis += "))";

					y_axis += std::to_string(j);
					y_axis += ")";

					y_value += std::to_string(Quals_ByPosition[i][j] / read_counter );
					y_value += ")";
					}
				else{
					label += name;
					label += ",";
					
					x_axis += std::to_string(i+1);
					x_axis += ",";

					y_axis += std::to_string(j);
					y_axis += ",";

					y_value += std::to_string(Quals_ByPosition[i][j] / read_counter );
					y_value += ",";
					}
				}
			else{
				label += name;
				label += ",";
				
				x_axis += std::to_string(i+1);
				x_axis += ",";

				y_axis += std::to_string(j);
				y_axis += ",";

				y_value += std::to_string(Quals_ByPosition[i][j] / read_counter );
				y_value += ",";		
				}
 		}
	}


R << "df3 <- data.frame (" << x_axis << " , "<< y_axis << " , "<< y_value << " , "<< label <<")"<< std::endl;
R <<  "ggplot(df3, aes(x= x_label,y=y_label,fill=labels, weight = y_weights)) ";
R << "+ geom_violin(scale=\"count\",position=\"identity\",linetype=\"blank\") ";
R << "+theme(axis.text.x = element_text(size=4,angle=90),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
R << "+ ggtitle(\"Quality Score Distributions by Position\")+ xlab(\"Position\") ";
R << "+ ylab(\"Quality Scores\") + labs(\""<< args.basename <<"\")" << std::endl;




/*passing reads dump*/
std::string value = "value = c(";
y_axis  = "y_label = factor(c(";
x_axis  = "x_label = factor(c(";

int passing_reads_size = Passing_Reads.size();

for (int i =0; i < passing_reads_size; ++i){
	for (int j = 0; j <= args.max_quality_score; ++j){
		if (Passing_Reads[i][j] == 0){
			value+= "NA";
			}
		else{
			value += std::to_string(Passing_Reads[i][j]/read_counter);
			}
		x_axis+=std::to_string(i+1);
		y_axis+=std::to_string(j);

		if (i == passing_reads_size-1){
			if (j == args.max_quality_score ){
				y_axis+="))";
				x_axis += "))";
				value+=")";
				}
			else{
				y_axis+=",";
				x_axis += ",";	
				value+=",";
				}
			}
		else{
			y_axis+=",";
			x_axis += ",";
			value+=",";
			}			
		}
	}


R << "df4 <- data.frame (" << x_axis << " , "<< y_axis << " , "<< value <<")"<< std::endl;
R << "ggplot(df4, aes(x = x_label, y = y_label)) + geom_tile(aes(fill = value)) ";
R << " + scale_fill_gradient(low=\"steelblue\", high=\"white\")+theme(axis.text.x = element_text(size=4,angle=90),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())";
R << " + ggtitle(\"Passing reads for "<<args.basename<<"\") ";
R << "  + ylab(\"Minimum Quality Score\")+xlab(\"Number of bases >= Minimum Quality\")";
R << " +labs(fill = \"Proportion of Reads\")" << std::endl;


return 0;
}
