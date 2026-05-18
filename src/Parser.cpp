#include "Parser.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>


namespace config{
    using namespace std;

    vector<Nucleotide> Parser::parse_sequence(string& file_path){
        ifstream file(file_path);
        if(!file.is_open()){
            throw runtime_error("Unable to open file");
        }

        vector<Nucleotide> sequence;
        string curr_line = "";

        while (getline(file, curr_line)){
            //Skip empty line
            if (curr_line.empty()){
                continue;
            }

            //Skip header line
            if (curr_line[0] == '>'){
                continue;
            }

            //Convert letters into Nucleotides
            for(char c : curr_line){
                c = toupper(c);
                switch(c){
                    case 'A': sequence.push_back(Nucleotide::A); break;
                    case 'C': sequence.push_back(Nucleotide::C); break;
                    case 'G': sequence.push_back(Nucleotide::G); break;
                    case 'T': sequence.push_back(Nucleotide::T); break;
                    default:
                        cout << "Invalid sequence - Nucleotide not A, C, G, or T" << endl;
                        break;
                }
            }
        }
        file.close();
        return sequence;
    }
}