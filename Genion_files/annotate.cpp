#include <utility>

#include <unordered_set>
#include <set>
#include <numeric>
#include <vector>
#include <tuple>
#include <functional>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include <cmath>

#include <cxxopts.hpp>
#include <IITree.h>

#include "annotate.h"
#include "candidate.h"
#include "locus.h"
#include "paf.h"
#include "util.h"
#include <stats/stattest.h>

#include <iostream>

#include <cassert>
#include <filesystem>

using std::string;
using std::ostream;
using std::vector;

template< class B>
inline void print_tsv(ostream &ost, B b){
    ost << b << "\n";
}

template <class B, class... A>
inline void print_tsv(ostream &ost, B b, A... a){
    ost << b << "\t";
    print_tsv(ost, a...);
}


namespace annotate{



    ostream& operator<<(ostream& os, const locus &lc){
        os << lc.chr << "\t" << lc.position;
        return os;
    }
    enum class SEQDIR{ forward, reverse, unknown};
    cxxopts::ParseResult parse_args(int argc, char **argv ){
        try{
            cxxopts::Options *options = new cxxopts::Options(argv[0], "Fusion candidate annotation");

            options->add_options()


                ("i,input", "Output path of Genion filter stage", cxxopts::value<string>())
                ("o,output", "Output path of Genion annotation stage", cxxopts::value<string>())
                ("s,minsupport", "min support to flag PASS", cxxopts::value<size_t>()->default_value("3"))
                ("maxrtfin", "maximum allowed fin for a read-through event, "
                                    "if larger event will be treated as an SV", cxxopts::value<double>()->default_value("0.5"))
                ("maxrtdistance", "maximum allowed distance for a read-through event, "
                                    "if larger event will be treated as an SV", cxxopts::value<long>()->default_value("600000"))
                ("d,duplications", "genomicSuperDups.txt, unzipped",cxxopts::value<string>())//can be found at http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz
                ("r,reference", "Reference path used in filter stage",cxxopts::value<string>())
                ("c,keep_non_coding", "Keep non coding genes", cxxopts::value<bool>()->default_value("false"))
                ("h,help", "Prints help")
                ;
            cxxopts::ParseResult result = options->parse(argc, argv);
            int ret = 0;
            if( result.count("h")){
             //   std::cerr << options->help() << std::endl;
                ret |=1;
            }

            if(!result.count("i")){
                std::cerr << "input is required" << std::endl;
                ret |=8;
            }
            if(!result.count("o")){
                std::cerr << "output is required" << std::endl;
                ret |=16;
            }

            if(!result.count("d")){
                std::cerr << "Duplication annotation is required" << std::endl;
                ret |=32;
            }
            if(!result.count("r")){
                std::cerr << "reference is required" << std::endl;
                ret |=64;
            }

            if(ret != 0 ){
                std::cerr  << "\n" << options->help() << std::endl;
                exit(-1);
            }
            return result;
        }
        catch (const cxxopts::OptionException& e)
        {
            std::cerr << "error parsing options: " << e.what() << std::endl;
            exit(1);
        }
    }

//585     chr1    10000   87112   chr15:101906152 0       -       chr15   101906152       101981189       75037   11764   1000    N/A     N/A     N/A     N/A     align_both/0009/both0046049     77880   71      3611  74269   73743   526     331     195     0.992918        0.991969        0.00711601      0.00711937 
    IITree<locus, std::tuple<string, int, int, double> > read_duplication_annotation(string path){
        IITree<locus, std::tuple<string, int, int, double> > duplications;

        std::ifstream dup_file(path);
        if(!dup_file.is_open()){
            std::cerr << "[ERROR] Cannot open file:" << path << std::endl;
            exit(-1);
        } 
        string line;
        
        string ch;
        int start;
        int end;

        string m_ch;
        int m_start;
        int m_end;

        double frac_match;
        while(std::getline(dup_file, line)){
            vector<string> fields = rsplit(line,"\t");
            ch = fields[1];
            start = stoi(fields[2]);
            end = stoi(fields[3]);

            m_ch = fields[7];
            m_start = stoi(fields[8]);
            m_end = stoi(fields[9]);

            frac_match = stof(fields[26]);

            if(ch.find("chr")!= string::npos){
                ch = ch.substr(3);
            }
            if(m_ch.find("chr")!= string::npos){
                m_ch = m_ch.substr(3);
            }
            locus s_s(ch,start);
            locus s_e(ch,end);
//            std::cerr << m_ch << "\t" << m_start << "\t" << m_end << "\t" << s_s.chr << "\t" << s_s.position << "\t" << s_e.chr << "\t" << s_e.position << "\tDUP" << "\n";
            duplications.add(s_s, s_e, std::make_tuple(m_ch,m_start,m_end,frac_match));
        }
        dup_file.close();
        duplications.index();
        return duplications;
    }
  
    class interval{

        public:
        string chr;
        int start;
        int end;
        bool reverse_strand;

        interval( const string &chr, int start, int end, bool reverse_strand) :
            chr(chr),
            start(start),
            end(end),
            reverse_strand(reverse_strand) {}
        std::pair<locus, locus> as_loci()const {
            return std::make_pair(locus(chr,start),locus(chr,end));
        }
        
        interval () : chr(""), start(-1), end(-1), reverse_strand(0) {}
        
        std::pair< locus,locus> as_loci(){
            return std::make_pair(locus(chr,start),locus(chr,end));
        }
        bool overlaps(const interval &other) const{
            if( chr != other.chr){
                return false;
            }
            if(start > other.start){
                return start < other.end;
            }
            return end > other.start;
        }

        bool operator==(const interval &other) const{
            return chr == other.chr && start == other.start && end == other.end && reverse_strand == other.reverse_strand;
        }

        bool extend(const interval& other){

            if( *this == interval{}){
                *this = other;
            }
            else{
                assert(chr == other.chr);
                
                if(reverse_strand != other.reverse_strand){
                    return false;        
                }
                start = std::min(start, other.start);
                end   = std::min(end, other.end);
            }
            return true;
        }
        friend ostream& operator<<(ostream& os, const interval &lc){
            os  << lc.chr << ":" << lc.start << "-" << lc.end << (lc.reverse_strand?"-":"+");
            return os;
        }
    };
//chr1    HAVANA          gene    65419   71585   .       +       .       gene_id "ENSG00000186092.6"; gene_type "protein_coding"; gene_name "OR4F5"; level 2; hgnc_id "HGNC:14825"; havana_gene "OTTHUMG00000001094.4";
//   1    ensembl_havana  gene    65419   71585   .       +       .       gene_id "ENSG00000186092"; gene_version "6"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding";`
    

    class gene{
        public:
        interval range;

        bool reverse_strand;
        string gene_id;
        string gene_name;
        string gene_type;
        bool coding;

        gene( interval range, const string &gene_id, const string &gene_name, 
                const string &gene_type) 
            : range(range),
              gene_id(gene_id),
              gene_name(gene_name),
              gene_type(gene_type),
              coding(gene_type=="protein_coding")
        {}
        gene( interval range, const string &gene_id, const string &gene_name, 
                const string &gene_type, bool coding) 
            : range(range),
              gene_id(gene_id),
              gene_name(gene_name),
              gene_type(gene_type),
              coding(coding)
        {}
        gene(){}
    };

    class exon{
        public:
        interval range;
        string gene_id;  
        string transcript_id;  
        int exon_no;

        exon( interval range, const string &gene_id, const string &transcript_id, int exon_no) :
            range(range),
            gene_id(gene_id),
            transcript_id(transcript_id),
            exon_no(exon_no) {}
        friend ostream& operator<<(ostream& os, const exon &lc){
            os <<  lc.gene_id << "\t"  << lc.transcript_id <<"\t"<< lc.exon_no <<"\t" << lc.range;
            return os;
        }
    };

    class candidate_read{
        public:

        string read_id;
        vector<std::pair<interval, exon> > blocks;
        candidate_read(const string &rid) : read_id(rid) {}
        candidate_read(const Candidate &cand) : read_id(cand.id){
            for(const auto &p: cand.canonical){
                add_block(p);
            }
        }
        vector<int> first_exons;

        auto ranges() const -> std::map<string, interval> {
            std::map<string, interval> gene_ranges;
            for(const auto &ie : blocks){
                const string &gene_id = ie.second.gene_id;
                gene_ranges[gene_id].extend(ie.first);
            }
            return gene_ranges;
        }
        void log(std::ofstream &ost) const{
            ost << read_id;
            for(const auto &p: ranges()){
                ost << "\t"<< p.first << "\t" << p.second.chr << ":" << p.second.start << "-" << p.second.end;
            }
            ost << "\n";
        }

        std::map<string, locus> get_breakpoints(bool direction) const {
            
            std::map<string, locus> bps;
            
            string first_gene = blocks[0].second.gene_id;
            for(auto &block : blocks){

                string gene_id = block.second.gene_id;
                bool is_first = (gene_id == first_gene) == direction;

                bool reverse = block.first.reverse_strand;
                auto loci = block.first.as_loci();
                
                auto gene_ptr  = bps.find(gene_id);
                if ( gene_ptr == bps.end()){
                    if(reverse != is_first){
                        bps.emplace(gene_id,loci.first);
                    }
                    else{
                        bps.emplace(gene_id,loci.second);
                    }
                }
                else{
                    int pos = gene_ptr->second.position;

                    if(reverse){
                        if(is_first){
                            if(pos < loci.second.position){
                                bps[gene_id] = loci.second;
                            }
                        }
                        else{
                            if(pos > loci.first.position){
                                bps[gene_id] = loci.first;
                            }
                        }
                    }
                    else{
                        if(is_first){
                            if(pos > loci.first.position){
                                bps[gene_id] = loci.first;
                            }

                        }
                        else{
                            if(pos < loci.second.position){
                                bps[gene_id] = loci.second;
                            }

                        }
                    }
                }
            }

            return bps;
        }
        void add_block(const string &line){
            vector<string> fields = rsplit(line, "\t");

            int start = stoi(fields[1]);
            int end   = stoi(fields[2]);
            string chr = fields[3];
            bool reverse_strand = fields[6] == "1";

            int ex_start = stoi(fields[8]);
            int ex_end   = stoi(fields[9]);
            
            bool ex_rev_strand = fields[10] == "1";
            string gene_id  = fields[11];
            string transcript_id = fields[12];
            int exon_no               = stoi(fields[13]);
            if(exon_no == 1){
                first_exons.push_back(blocks.size());
            }
            
            
            interval alig(chr,start,end, reverse_strand);

            interval expos(chr,ex_start,ex_end, ex_rev_strand);
            exon ex(expos, gene_id, transcript_id, exon_no);

            blocks.push_back(std::make_pair(alig,ex));
        }
        void add_block(const std::pair<aligned_segment,::exon> &p){

            interval alig(p.first.chr, p.first.tmplt.start, p.first.tmplt.end, p.first.reverse_complemented);

            interval expos(p.second.chr, p.second.start, p.second.end, p.second.strand);
            exon ex(expos, p.second.gene_id, p.second.transcript_id, p.second.exon_number);
            if(p.second.exon_number == 1){
                first_exons.push_back(blocks.size());
            }

            blocks.push_back(std::make_pair(alig,ex));
        }
    };


    auto median(const vector<int> &values) -> double{
        int i = values.size() / 2;
        if(values.size() % 2 == 0){
            return values[i] / 2.0 + values[i + 1] / 2.0;
        }
        else{
            return values[i + 1];
        }
    }

    class candidate_fusion{

        public:
        std::map<string,double> non_covered_sum_ratio;

        string name;
        string id;
        vector<candidate_read> forward;
        vector<candidate_read> backward;

        vector<candidate_read> no_first;
        vector<candidate_read> multi_first;

        vector<std::pair<interval, interval> > duplications;
        vector<std::pair<gene, gene> > gene_overlaps;
        int invalid {0};


        void log(std::ofstream &ost) const {

            for(const candidate_read &cr : forward){
                cr.log(ost);
            }
            for(const candidate_read &cr : backward){
                cr.log(ost);
            }
            for(const candidate_read &cr : no_first){
                cr.log(ost);
            }
            for(const candidate_read &cr : multi_first){
                cr.log(ost);
            }
        }

        auto median_range() const -> vector<std::tuple<string, int, int>>{
            vector<std::tuple<string, int, int>> median_values;
            std::map<string, vector< int>> begins;
            std::map<string, vector< int>> ends;
            std::map<string, string> chrs;
            for(const candidate_read &cr : forward){
                for( const auto &pp: cr.ranges()){
                    begins[pp.first].push_back(pp.second.start);
                    ends[pp.first].push_back(pp.second.end);
                    chrs[pp.first] = pp.second.chr;
                }
            }
            for(const candidate_read &cr : backward){
                for( const auto &pp: cr.ranges()){
                    begins[pp.first].push_back(pp.second.start);
                    ends[pp.first].push_back(pp.second.end);
                    chrs[pp.first] = pp.second.chr;
                }
            }
            for(const candidate_read &cr : no_first){
                for( const auto &pp: cr.ranges()){
                    begins[pp.first].push_back(pp.second.start);
                    ends[pp.first].push_back(pp.second.end);
                    chrs[pp.first] = pp.second.chr;
                }
            }
            for(const candidate_read &cr : multi_first){
                for( const auto &pp: cr.ranges()){
                    begins[pp.first].push_back(pp.second.start);
                    ends[pp.first].push_back(pp.second.end);
                    chrs[pp.first] = pp.second.chr;
                }
            }

            for( auto &pp: begins){
                const auto &gn = pp.first;
                auto &bvec = pp.second;
                auto &evec = ends[gn];
                const string &chr = chrs[gn];
                sort(bvec.begin(), bvec.end());
                sort(evec.begin(), evec.end());
                median_values.emplace_back(chr, median(bvec), median(evec));
            }
            return median_values;
        }
        size_t total_count() const {
            return forward.size() 
                + backward.size() 
                + multi_first.size() 
                + no_first.size();
        }
        std::map<string, interval> fusion_gene_intervals(){
            std::map<string, string> chrs;
            std::map<string, int> mins; 
            std::map<string, int> maxs;
            std::map<string, bool> rev;
            for(const auto &v : {forward, backward, no_first, multi_first}){
                for(const auto &c : v){
                    for(const auto &i_e : c.blocks){
                        const interval &i = i_e.first;
                        const exon &e = i_e.second;
                        int mn = mins[e.gene_id];
                        int mx = maxs[e.gene_id];
                        if(i.start < mn || mn == 0){
                            mins[e.gene_id] = i.start;
                        }
                        if(i.end > mx){
                            maxs[e.gene_id] = i.end;
                        }
                        chrs[e.gene_id] = i.chr;
                        rev[e.gene_id] = i.reverse_strand;
                    }
                }
            } 
            std::map<string, interval> ivals;
            for(const auto &k_v : mins){
                const string &key = k_v.first;
                const int &mn = k_v.second;
                const int &mx = maxs[key];
                const string &chr = chrs[key];
                bool rs = rev[key];
                ivals.emplace(key, interval{chr,mn,mx,rs});
            }
            return ivals;
        }


        candidate_fusion() {}
    };

    auto dash_fold(const string &a, const string &b){
        return std::move(a) + "::" + b;
    }
    class fusion_manager{
        public:
        std::map<string, candidate_fusion> fusions;
        std::map<string, int> gene_counts;

        fusion_manager( const vector<Candidate> &candidates, const std::unordered_map<string, gene> &gene_annot,
            const std::unordered_map<string, int> &exon_counts){
            for( auto &cand : candidates){

                add_read(candidate_read{cand}, gene_annot, exon_counts);
            }

        }
        fusion_manager() {}
        void add_read(const candidate_read &read, const std::unordered_map<string, gene> &gene_annot,
            const std::unordered_map<string, int> &exon_counts){


            std::set<string> gene_ids;
            std::map<string,int> gene_order;   //use
            std::map<string, std::unordered_set<string>> transcript_ids;
            std::map<string, double> approximate_coverage;
            int and_all_blocks  = 1;
            int not_and_all_blocks = 1;
            int index = 0;
            for(auto i_and_e : read.blocks){
                auto ite = gene_ids.find(i_and_e.second.gene_id);
                if(ite == gene_ids.end()){
                    gene_order[i_and_e.second.gene_id] = index;
                    ++index;
                }
                gene_ids.insert(i_and_e.second.gene_id);
                
                bool exon_strand = i_and_e.second.range.reverse_strand;
                bool interval_strand = i_and_e.first.reverse_strand;

                int strand_xor = exon_strand ^ interval_strand;

                and_all_blocks = and_all_blocks && strand_xor;
                not_and_all_blocks = not_and_all_blocks && (! strand_xor);

                
                if(gene_annot.find(i_and_e.second.gene_id) == gene_annot.end()){
                    std::cerr << i_and_e.second.gene_id << " is not in annotation!\n";
                }
                //int exon_count = exon_counts.at(i_and_e.second.transcript_id);
                transcript_ids[i_and_e.second.gene_id].insert(i_and_e.second.transcript_id);
                approximate_coverage[i_and_e.second.gene_id] += 1;//(1.0/exon_count);
            }
            

            string fusion_name = "";
            for(const string &id : gene_ids){
                fusion_name += gene_annot.find(id)->second.gene_name + "::";
            }

            fusion_name.pop_back();
            fusion_name.pop_back();

            string fusion_id = std::accumulate( std::next(std::begin(gene_ids)), std::end(gene_ids), *(std::begin(gene_ids)), dash_fold);

            for( string gid : gene_ids){
                gene_counts[gid]+=1;
            }
            auto &cand = fusions[fusion_id];
            if( !(and_all_blocks || not_and_all_blocks)){
                //Invalid Strand configuration
                //std::cerr << read.read_id << "\n";
                //for(const string &id : gene_ids){
                //    std::cerr << id<< "::";
               // }
                cand.invalid +=1;
            }

            for( string gid : gene_ids){
                int max_exon_count = 1;
                for(string tid : transcript_ids[gid]){
                    int exon_count = exon_counts.at(tid);

                    if( max_exon_count < exon_count){
                        max_exon_count = exon_count;
                    }
                }
                //std::cerr << gid << "\t" << max_exon_count << "\t" << approximate_coverage[gid] << "\n";
                cand.non_covered_sum_ratio[gid]+= 10.0 / (10 + max_exon_count - approximate_coverage[gid]);
            }
                
            cand.name = fusion_name;
            cand.id = fusion_id;
            int last_first = - 1;
            if(read.first_exons.size() > 1){
                cand.multi_first.push_back(read);
                return;
            }
            if(read.first_exons.size() == 0){
                cand.no_first.push_back(read);
                return;
            }
            last_first = read.first_exons.back();
            if(read.blocks[last_first].second.gene_id == *(gene_ids.rbegin())){
                cand.forward.push_back(read);
            }
            else{
                cand.backward.push_back(read);
            }
        }
    };

    bool make_gene(const string &line, gene &g){
        vector<string> tabs = rsplit(line, "\t");
        string ch(tabs[0]);
        if(ch.find("chr")!=string::npos){
            ch = ch.substr(3);
        }
        interval range( tabs[0], stoi(tabs[3]), stoi(tabs[4]), tabs[6]=="-");
        if(tabs[2] != "gene"){
            return false;
        }
        vector<string> fields = rsplit(tabs[8], ";");
        string gene_id, gene_name, gene_type;

        for(auto iter = fields.begin(); iter != fields.end(); iter++){

            if( iter->find("gene_id") != string::npos){
                string _id = iter->substr(iter->find("d ")+3);
                _id.pop_back();
                if( _id.find(".") != string::npos){
                    size_t dot_pos = _id.find(".");
                    _id = _id.substr(0,dot_pos);
                }
                gene_id = _id;
            } 

            if( iter->find("gene_name") != string::npos){
                string _id = iter->substr(iter->find("e ")+3);
                _id.pop_back();
                gene_name = _id;
            }
            if( iter->find("gene_biotype") != string::npos ||
                iter->find("gene_type") != string::npos 
              ){
                string _id = iter->substr(iter->find("e ")+3);
                _id.pop_back();
                gene_type = _id;
            }
        }
   //     std::cerr << gene_type << "\tTYPE\n";
        g = gene(range, gene_id, gene_name, gene_type);
        return true;
    }

// umap of transcript to last exon id
    std::unordered_map<string, int> read_last_exons(string gtf_path){
        std::ifstream gtf_file(gtf_path);
        if(!gtf_file.is_open()){
            std::cerr << "[ERROR] Cannot open file:" << gtf_path << std::endl;
            exit(-1);
        } 


        std::unordered_map<string, int> int_map;
        string line;
        while(std::getline(gtf_file, line)){
            if(line[0]=='#'){ //Comment
                continue;
            }
            vector<string> tabs = rsplit(line, "\t");
            if(tabs[2] != "exon"){
                continue;
            }

            vector<string> fields = rsplit(tabs[8], ";");
            string transcript_id = "-1";
            int exon_number = -1;
            for(auto iter = fields.begin(); iter != fields.end(); iter++){

                if( iter->find("transcript_id") != string::npos){
                    string _id = iter->substr(iter->find("d ")+3);
                    _id.pop_back();
                    if( _id.find(".") != string::npos){
                        size_t dot_pos = _id.find(".");
                        _id = _id.substr(0,dot_pos);
                    }
                    transcript_id = _id;

                }
                if( iter->find("exon_number") != string::npos){
                    string _id = iter->substr(iter->find("r ")+3);
                    _id.pop_back();
                    if( _id.find(".") != string::npos){
                        size_t dot_pos = _id.find(".");
                        _id = _id.substr(0,dot_pos);
                    }
                    exon_number = stoi(_id);
                        
                }
            }
            if( transcript_id == "-1"){
                std::cerr << "Transcript doesn't have transcript_id\n";
            }
            if( exon_number > int_map[transcript_id]){
                int_map[transcript_id] = exon_number;
            }
        }
        gtf_file.close();

        return int_map;
    }

    template<class K, class V>
    vector<std::pair<K, K>> get_key_pairs( const std::map<K,V> &map){
        vector<std::pair<K,K>> pairs;
        for(auto iter = std::begin(map); iter != std::end(map); ++iter){
            for(auto inner = std::next(iter); inner !=std::end(map); ++inner){
                pairs.emplace_back(iter->first, inner->first);
            }
        }
        return pairs;
    }

    std::unordered_map<string, int> read_transcript_exon_counts(string gtf_path){
        std::ifstream gtf_file(gtf_path);
        if(!gtf_file.is_open()){
            std::cerr << "[ERROR] Cannot open file:" << gtf_path << std::endl;
            exit(-1);
        } 

        std::unordered_map<string, int> transcript_exon_counts;
        string line;
        while(std::getline(gtf_file, line)){
            if(line[0]=='#'){ //Comment
                continue;
            }

            vector<string> tabs = rsplit(line, "\t");
            string ch(tabs[0]);
            if(ch.find("chr")!=string::npos){
                ch = ch.substr(3);
            }
            interval range( tabs[0], stoi(tabs[3]), stoi(tabs[4]), tabs[6]=="-");
            if(tabs[2] != "exon"){
                continue;
            }
            vector<string> fields = rsplit(tabs[8], ";");
            string transcript_id;
            for(auto iter = fields.begin(); iter != fields.end(); iter++){


                if( iter->find("transcript_id") != string::npos){
                    string _id = iter->substr(iter->find("d ")+3);
                    _id.pop_back();
                    if( _id.find(".") != string::npos){
                        size_t dot_pos = _id.find(".");
                        _id = _id.substr(0,dot_pos);
                    }
                    transcript_id = _id;
                    transcript_exon_counts[transcript_id]+=1;
                    break;
                } 
/*  In case we have to handle duplicate exons.
                if( iter->find("exon_number") != string::npos){
                    string _id = iter->substr(iter->find("d ")+3);
                    _id.pop_back();
                    if( _id.find(".") != string::npos){
                        size_t dot_pos = _id.find(".");
                        _id = _id.substr(0,dot_pos);
                    }
                    exon_number = _id;
*/
                } 
            }
            
        gtf_file.close();
        return transcript_exon_counts;
    }

    std::unordered_map<string, gene> read_gene_annotation(string gtf_path){
        std::ifstream gtf_file(gtf_path);
        if(!gtf_file.is_open()){
            std::cerr << "[ERROR] Cannot open file:" << gtf_path << std::endl;
            exit(-1);
        } 

        std::unordered_map<string, gene> valid_set;
        string line;
        while(std::getline(gtf_file, line)){
            if(line[0]=='#'){ //Comment
                continue;
            }
            gene g;
            if(make_gene(line, g)){
                valid_set.emplace( g.gene_id, g);
            }
        }
        gtf_file.close();
        return valid_set;
    }

    void annotate_duplications_and_overlaps(fusion_manager &fm,
            const std::unordered_map<string, gene> &gene_annot,
            const string &dup_path){

        IITree<locus, std::tuple<string, int, int, double> > duplications = read_duplication_annotation(dup_path);;
        vector< size_t> overlaps;
        for( auto &cand : fm.fusions){
            std::map<string, interval> ivals = cand.second.fusion_gene_intervals();
            auto key_pairs = get_key_pairs(ivals);
                
            //Duplication annotation
            for(const auto &key_pair : key_pairs){
                const string &f = key_pair.first;
                const string &s = key_pair.second;
                
                const interval &i = (ivals.find(f))->second;
                const auto loci = i.as_loci();

                duplications.overlap(loci.first,loci.second, overlaps);

                for(size_t d : overlaps){
                    const auto &dup = duplications.data(d);
                    interval l(std::get<0>(dup), std::get<1>(dup), std::get<2>(dup), 0);
                    interval &r =( ivals.find(s))->second; 
                    if(l.overlaps(r)){
                        cand.second.duplications.emplace_back(i,l);
                    }
                }
                overlaps.clear(); 
            }
            //X

            //Gene overlap annotation
            for(const auto &key_pair : key_pairs){
                const gene &f = gene_annot.find(key_pair.first)->second;
                const gene &s = gene_annot.find(key_pair.second)->second;
                if(f.range.overlaps(s.range)){
                    cand.second.gene_overlaps.emplace_back(f,s); 
                }             
            }
            //X
        }
    }


    std::pair<size_t,size_t> count_genes( const string &feature_table_path,
            std::unordered_map<string, size_t> &count_table,
            bool all = false){
        
          
        
        std::ifstream feature_file(feature_table_path);
        
        size_t total_normal_count = 0;
        size_t total_chimer_count = 0;

        string line;
        while(std::getline(feature_file, line)){
            vector<string> fields =  rsplit(line, "\t");
            if(!all && stoi(fields[2]) !=0){ // Split Alignment
                ++total_chimer_count;
                continue;
            }
            ++total_normal_count;
            
            string gene_id1 = fields[1].substr(0,15);
            string gene_id2 = fields[1].substr(17);
            if(all){
                count_table[gene_id1]+=1;
                if(gene_id1 != gene_id2){
                    count_table[gene_id2]+=1;
                }
            }
            else{
                if(gene_id1 == gene_id2){
                    count_table[gene_id1]+=1;
                }
            }

        }
        return std::make_pair(total_normal_count,total_chimer_count);
    }


    std::unordered_map<string, SEQDIR>  read_read_directions(const string &path){
        
        std::unordered_map<string, SEQDIR> directions; 
        std::ifstream dir_file(path);
        string line;

        while(std::getline(dir_file, line)){
            vector<string> fields = rsplit(line,"\t");
            if(fields[1] =="NONE"){
                directions[fields[0]] = SEQDIR::unknown;
                continue;
            }
            if(fields[1] == "A" && stoi(fields[2]) > 50){
                directions[fields[0]] = SEQDIR::reverse;
            }
            else if (fields[1] == "T" && stoi(fields[2]) < 50){
                directions[fields[0]] = SEQDIR::forward;
            }
            else{
                directions[fields[0]] = SEQDIR::unknown;
            }
        }
        return directions;
    }

    bool is_cluster_rt( const candidate_fusion &cf, double fin,
            double forw_rt_ex, double back_rt_ex,
            int max_rt_distance = 50000, double max_fin = 0.1){
        
        const candidate_read *read = NULL;
        if( cf.forward.size() > 0){
            read = &cf.forward[0];
        }
        else if( cf.backward.size()){
            read = &cf.backward[0];
        }
        else if( cf.multi_first.size()){
            read = &cf.multi_first[0];
        }
        else if( cf.no_first.size()){
            read = &cf.no_first[0];
        }
        else{
            return false;
        }
        //vector<std::pair<interval, exon> > blocks;
    
        if( read->blocks.size() < 2){
            return false;
        }
        auto &blocks = read->blocks;
        size_t i = 0;
        for( i=1; i < blocks.size(); ++i){
            if( blocks[i].second.gene_id != blocks[i-1].second.gene_id){
                break;
            }
        }

        auto b1 = blocks[i-1];
        auto b2 = blocks[i];

        if( b1.second.range.chr != b2.second.range.chr){
            return false;
        }
        std::array<int, 4> positions {{ b1.second.range.start, b1.second.range.end,  b2.second.range.start, b2.second.range.end}};
        std::sort(positions.begin(),positions.end());
        int distance = positions[2] - positions[1];
        if( distance > max_rt_distance){
            return false;
        }
        if( forw_rt_ex < 0.8 || back_rt_ex < 0.8){
            return false;
        }
        if( fin > max_fin){
            return false;
        }

        return true;
    }
   
    double statistically_test_candidate(const candidate_fusion &fusion,
            double chimera_rate,
            const std::unordered_map<string, size_t> &gene_counts
            ){

            
        vector<string> genes = rsplit(fusion.id,"::");
        vector<size_t> normal_counts;
        for(const string &gene : genes){

            auto count_ptr = gene_counts.find(gene);
            size_t count = 0;
            if(count_ptr != gene_counts.end()){
                count = count_ptr->second;
            }
            normal_counts.push_back(count);
        }

        size_t mult_count = std::accumulate(normal_counts.begin(), normal_counts.end(), 1L, std::multiplies<size_t>());
        double average_normal_count = std::pow(static_cast<double>(mult_count), 1.0/normal_counts.size());

        int x = fusion.total_count();
        int n = x + average_normal_count;
        int m = x + chimera_rate * average_normal_count;
        int N = 2 * n; 
        double pvalue = hyper_geom_cdf(x,n,m,N);
        
        return pvalue;
    }

    int annotate_calls_direct( 
            const string &output_path,
            const string &log_path,
            const string &gtf_path,
            const string &duplication_path,
            const vector<Candidate> &candidates,
            std::unordered_map<string, size_t> gene_counts,
            size_t min_support,
            int total_normal_count, int total_chimer_count,
            int maxrtdistance,      double maxrtfin,
            bool only_coding){

        bool full_debug_output = false;

        std::unordered_map<string, gene> gene_annot = read_gene_annotation(gtf_path);
        
        std::unordered_map<string, int> last_exons = read_last_exons(gtf_path);        //May be removed.
        std::unordered_map<string, int> transcript_exon_counts = read_transcript_exon_counts(gtf_path);
        /*
        vector<candidate_read> candidate_reads;
        for( const Candidate &cand: candidates){
            candidate_read cr(cand.id);
            for(const auto &p: cand.canonical){
                cr.add_block(p);
            }
            candidate_reads.push_back(cr);
        };
       */ 
        fusion_manager fm{candidates, gene_annot, transcript_exon_counts};
//        for( auto &cand : candidate_reads){
//            fm.add_read(cand, gene_annot, transcript_exon_counts);
//        }
        
        annotate_duplications_and_overlaps(fm, gene_annot, duplication_path);
        

        double mean_chimera_ratio = static_cast<double>(total_chimer_count) / total_normal_count;


        vector<double> pvalues;
        for( const auto &cand : fm.fusions){
            double pvalue = statistically_test_candidate(cand.second, mean_chimera_ratio, gene_counts);
            pvalues.push_back(pvalue);
        }
        auto hypothesis = multiple_test(pvalues, 0.05, pvalue_corrector::BENJAMINI_YEKUTIELI);



        std::ofstream outfile( output_path );
        std::ofstream logfile( log_path );
        std::ofstream outfile_fail( output_path + ".fail");

        auto pval_iter = pvalues.begin();
        auto corr_pval_iter = hypothesis.corr_pvals.begin();
        auto null_iter = hypothesis.null_rejected.begin();
        for( const auto &cand : fm.fusions){
            string fusion_id = cand.first;
            bool null_rejected = *null_iter;
            double pvalue = *pval_iter;
            double corr_pvalue = *corr_pval_iter;
            pval_iter = std::next(pval_iter);
            corr_pval_iter = std::next(corr_pval_iter);
            null_iter = std::next(null_iter);
            int total_count = cand.second.total_count();

            int total_count_putative_full_length = cand.second.forward.size() 
                + cand.second.backward.size();
            vector<string> genes = rsplit(fusion_id,"::");

            bool coding_flag = false;
            if( only_coding){
                for( const string &g : genes){
                    auto gptr = gene_annot.find(g);
                    if(gptr == gene_annot.end()){
                        std::cerr << "Gene " << g << " is not in annotation!\n";
                        continue;
                    }
                    if(gptr->second.coding == false){
                        coding_flag = true;
                        break;
                    }
                }
            }
            double gene_count_sum = 0;
            string gene_count_string = "";
            string idf_string = "";
            double total_idf = 0;
            for(const string &gene : genes){
                gene_count_sum += gene_counts[gene];
                gene_count_string+= std::to_string(gene_counts[gene]) + ";";
                idf_string+= std::to_string(fm.gene_counts[gene]-total_count) + ";";
                total_idf+= fm.gene_counts[gene] - total_count;
            }
            
            double tfidf_score = total_count * std::log(fm.fusions.size()/(1+total_idf/2));
            double tfidf_score_full_len = total_count_putative_full_length * std::log(fm.fusions.size()/(1+total_idf/2));
            
            int tcpflnz;
            if(total_count == 0){
                tcpflnz = 1;
            }
            else{
                tcpflnz = total_count;
            }

            double fin_score = genes.size() * total_count / (gene_count_sum+1);

            double fg_count = cand.second.non_covered_sum_ratio.at(genes[0]);
            double lg_count = cand.second.non_covered_sum_ratio.at(genes[1]);
            double forward_rt_ex  =  1.0 * fg_count / tcpflnz;
            double backward_rt_ex = 1.0 * lg_count / tcpflnz;
            double bad_strand_ratio = static_cast<double>(cand.second.invalid)/cand.second.total_count();
            string pass_fail_code = "";
            if(coding_flag){
                pass_fail_code += ":noncoding";
            }

            if(cand.second.gene_overlaps.size() > 0){
                pass_fail_code += ":overlaps";
            }
            if(cand.second.duplications.size() > 0){
                pass_fail_code += ":segdup";
            }
            if(bad_strand_ratio > 0.25){
                pass_fail_code += ":badstrand";
            }
            if( cand.second.forward.size() + cand.second.backward.size() 
                    + cand.second.multi_first.size() < min_support){
                pass_fail_code += ":lowsup";
            }
            if( pass_fail_code != ""){
                pass_fail_code = "FAIL" + pass_fail_code;
            }
            else{
                if( is_cluster_rt( cand.second, fin_score, forward_rt_ex, backward_rt_ex, 
                            maxrtdistance, maxrtfin)){
                    pass_fail_code = "PASS:RT";
                }
                else if( null_rejected){
                    pass_fail_code = "PASS:GF";
                }
                else{
                    pass_fail_code = "FAIL:RP";
                }
            }

/*
            std::map<string, vector<locus>> breakpoints;

            bool is_forward = true;
            for(const auto &ff : {cand.second.forward, cand.second.backward, cand.second.no_first, cand.second.multi_first}){
                for(const auto &fus : ff){
                    for(const auto &bp_pair : fus.get_breakpoints(is_forward)){
                        breakpoints[bp_pair.first].push_back(bp_pair.second);
                    }
                }
                is_forward = false;
            }

            std::map<string, std::pair<double, double>> breakpoint_ranges;
            std::map<string, string> chromosome_per_gene;
            for( const auto &bpp : breakpoints){
                chromosome_per_gene[bpp.first] = bpp.second[0].chr;
                breakpoint_ranges[bpp.first] = mean_and_std(bpp.second, [] (const locus &l) -> double {
                    return static_cast<double>(l.position);
                });
            }
            */
            if(full_debug_output){ 
            //#FusionID(Ensembl) Forward-Support Backward-Support Multi-First-Exon No-First-Exon Genes-Overlap Segmental-Duplication-Count FusionName(Symbol) FiN-Score Pass-Fail-Status total-normal-count fusion-count normal-counts proper-normal-count proper-FiN-Score total-other-fusion-count other-fusion-counts ffigf-score proper-ffigf-score A B Anorm Bnorm 
                outfile << fusion_id << "\t" << cand.second.forward.size() << "\t"
                    << cand.second.backward.size()  << "\t"
                    << cand.second.multi_first.size() << "\t" << cand.second.no_first.size()
                    << "\t" <<  cand.second.gene_overlaps.size() 
                    << "\t" <<  cand.second.duplications.size()
                    << "\t" << cand.second.name << "\t" << fin_score
                    << "\t" <<  pass_fail_code
                    << "\t" << gene_count_sum << "\t" << total_count <<  "\t"  << gene_count_string << "\t"
                    << total_count_putative_full_length << "\t" << genes.size() * total_count_putative_full_length / ( gene_count_sum + 1)
                    << "\t" <<  total_idf << "\t" << idf_string << "\t" << tfidf_score << "\t" << tfidf_score_full_len
                    << "\t" << fg_count  << "\t" << lg_count << "\t"
                    << forward_rt_ex << "\t" << backward_rt_ex << "\t"
                    << pvalue << "\t" << corr_pvalue << "\t" << (null_rejected?"pPASS":"pFAIL") 
                    << "\t" << static_cast<double>(cand.second.invalid)/cand.second.total_count() << "\n";
            }
            else{
                if(pass_fail_code.find("PASS")!=string::npos){
/*
                    string pos_string;
                    for( const string &g: genes){
                        auto rang = breakpoint_ranges.at(g);
                        string ch = chromosome_per_gene.at(g);
                        pos_string += g + "(" + ch + ":" + std::to_string(rang.first) + "±" + std::to_string(rang.second) + ")";
                    }

                    */
                    std::stringstream range_stream;
                    for(const auto &tup :cand.second.median_range()){
                        range_stream << std::get<0>(tup) << ":" << std::get<1>(tup) << "-" << std::get<2>(tup) << ";";
                    }
                    print_tsv(outfile, fusion_id, cand.second.name, tfidf_score_full_len, fin_score, total_count, gene_count_string, pass_fail_code, range_stream.str());// pos_string);

                    cand.second.log(logfile);
                }
                else{
                    print_tsv(outfile_fail, fusion_id, cand.second.name, tfidf_score_full_len, fin_score, total_count, gene_count_string, pass_fail_code);

                }
            }
        }
       
        logfile.close();
        outfile_fail.close();
        // Remove empty files if they exist and are empty
        if (std::filesystem::exists(log_path) && std::filesystem::file_size(log_path) == 0) {
            std::filesystem::remove(log_path);
        }
        if (std::filesystem::exists(output_path + ".fail") && std::filesystem::file_size(output_path + ".fail") == 0) {
            std::filesystem::remove(output_path + ".fail");
        }
        return 0;  
    }
       

    int annotate_calls(int argc, char **argv){
        auto  opt = parse_args(argc, argv);

        size_t min_support = opt["minsupport"].as<size_t>();

        string input_prefix(opt["input"].as<string>());
        string reference_path(opt["reference"].as<string>());

        bool filter_non_coding = !opt["c"].as<bool>();

        string gtf_path = reference_path + "/1.gtf";
        std::unordered_map<string, gene> gene_annot = read_gene_annotation(gtf_path);
        

        std::unordered_map<string, int> last_exons = read_last_exons(gtf_path);        //May be removed.
        std::unordered_map<string, int> transcript_exon_counts = read_transcript_exon_counts(gtf_path);

        string chains_path = input_prefix + "/chains.fixed.txt";

        vector<candidate_read> candidates;

        std::ifstream chain_file(chains_path);
        
        string line;
        while(std::getline(chain_file, line)){
            vector<string> fields =  rsplit(line, "\t");
            int block_count = stoi(fields[1]);
            candidate_read cr(fields[0]);
            for(int i = 0; i < block_count; ++i){
                std::getline(chain_file, line);
                cr.add_block(line);
            }
            candidates.push_back(cr);
        }

        chain_file.close();
        fusion_manager fm;
        for( auto &cand : candidates){
            fm.add_read(cand, gene_annot, transcript_exon_counts);
            //fm.add_read(cand, gene_annot, last_exons, read_directions);
        }
        
        annotate_duplications_and_overlaps(fm, gene_annot, opt["duplications"].as<string>());
        
        string feature_table_path = input_prefix + "/feature_table.tsv";

        std::unordered_map<string, size_t> gene_counts;
        auto[total_normal_count,total_chimer_count] = count_genes(feature_table_path, gene_counts, false);
        double mean_chimera_ratio = static_cast<double>(total_chimer_count) / total_normal_count;


        vector<double> pvalues;
        for( const auto &cand : fm.fusions){
            double pvalue = statistically_test_candidate(cand.second, mean_chimera_ratio, gene_counts);
            pvalues.push_back(pvalue);
        }
        auto hypothesis = multiple_test(pvalues, 0.05, pvalue_corrector::BENJAMINI_YEKUTIELI);

        auto pval_iter = pvalues.begin();
        auto corr_pval_iter = hypothesis.corr_pvals.begin();
        auto null_iter = hypothesis.null_rejected.begin();
        for( const auto &cand : fm.fusions){
            bool null_rejected = *null_iter;
            double pvalue = *pval_iter;
            double corr_pvalue = *corr_pval_iter;
            pval_iter = std::next(pval_iter);
            corr_pval_iter = std::next(corr_pval_iter);
            null_iter = std::next(null_iter);
            int total_count = cand.second.total_count();

            int total_count_putative_full_length = cand.second.forward.size() 
                + cand.second.backward.size();
            vector<string> genes = rsplit(cand.first,"::");

            bool coding_flag = false;
            if( filter_non_coding){
                for( const string &g : genes){
                    auto gptr = gene_annot.find(g);
                    if(gptr == gene_annot.end()){
                        std::cerr << "Gene " << g << " is not in annotation!\n";
                        continue;
                    }
                    if(gptr->second.coding == false){
                        coding_flag = true;
                        break;
                    }
                }
            }
            double gene_count_sum = 0;
            string gene_count_string = "";
            string idf_string = "";
            double total_idf = 0;
            for(const string &gene : genes){
                gene_count_sum += gene_counts[gene];
                gene_count_string+= std::to_string(gene_counts[gene]) + ";";
                idf_string+= std::to_string(fm.gene_counts[gene]-total_count) + ";";
                total_idf+= fm.gene_counts[gene] - total_count;
            }
            
            double tfidf_score = total_count * std::log(fm.fusions.size()/(1+total_idf/2));
            double tfidf_score_full_len = total_count_putative_full_length * std::log(fm.fusions.size()/(1+total_idf/2));
            
            int tcpflnz;
            if(total_count == 0){
                tcpflnz = 1;
            }
            else{
                tcpflnz = total_count;
            }

            double fin_score = genes.size() * total_count / (gene_count_sum+1);

            double fg_count = cand.second.non_covered_sum_ratio.at(genes[0]);
            double lg_count = cand.second.non_covered_sum_ratio.at(genes[1]);
            double forward_rt_ex  =  1.0 * fg_count / tcpflnz;
            double backward_rt_ex = 1.0 * lg_count / tcpflnz;
            double bad_strand_ratio = static_cast<double>(cand.second.invalid)/cand.second.total_count();
            string pass_fail_code = "";
            if(coding_flag){
                pass_fail_code += ":noncoding";
            }

            if(cand.second.gene_overlaps.size() > 0){
                pass_fail_code += ":overlaps";
            }
            if(cand.second.duplications.size() > 0){
                pass_fail_code += ":segdup";
            }
            if(bad_strand_ratio > 0.25){
                pass_fail_code += ":badstrand";
            }
            if( cand.second.forward.size() + cand.second.backward.size() 
                    + cand.second.multi_first.size() < min_support){
                pass_fail_code += ":lowsup";
            }
            if( pass_fail_code != ""){
                pass_fail_code = "FAIL" + pass_fail_code;
            }
            else{
                if( is_cluster_rt( cand.second, fin_score, forward_rt_ex, backward_rt_ex, 
                            opt["maxrtdistance"].as<long>(), opt["maxrtfin"].as<double>())){
                    pass_fail_code = "PASS:RT";
                }
                else if( null_rejected){
                    pass_fail_code = "PASS:GF";
                }
                else{
                    pass_fail_code = "FAIL:RP";
                }
            }
            
            //#FusionID(Ensembl) Forward-Support Backward-Support Multi-First-Exon No-First-Exon Genes-Overlap Segmental-Duplication-Count FusionName(Symbol) FiN-Score Pass-Fail-Status total-normal-count fusion-count normal-counts proper-normal-count proper-FiN-Score total-other-fusion-count other-fusion-counts ffigf-score proper-ffigf-score A B Anorm Bnorm 
            std::cout << cand.first << "\t" << cand.second.forward.size() << "\t"
                << cand.second.backward.size()  << "\t"
                << cand.second.multi_first.size() << "\t" << cand.second.no_first.size()
                << "\t" <<  cand.second.gene_overlaps.size() 
                << "\t" <<  cand.second.duplications.size()
                << "\t" << cand.second.name << "\t" << fin_score
                << "\t" <<  pass_fail_code
                << "\t" << gene_count_sum << "\t" << total_count <<  "\t"  << gene_count_string << "\t"
                << total_count_putative_full_length << "\t" << genes.size() * total_count_putative_full_length / ( gene_count_sum + 1)
                << "\t" <<  total_idf << "\t" << idf_string << "\t" << tfidf_score << "\t" << tfidf_score_full_len
                << "\t" << fg_count  << "\t" << lg_count << "\t"
                << forward_rt_ex << "\t" << backward_rt_ex << "\t"
                << pvalue << "\t" << corr_pvalue << "\t" << (null_rejected?"pPASS":"pFAIL") 
                << "\t" << static_cast<double>(cand.second.invalid)/cand.second.total_count() << "\n";

        }
       
        string bp_file_path = opt["output"].as<string>() + "/breakpoints.tsv";

        std::ofstream bp_file(bp_file_path);


        for( const auto &cand : fm.fusions){
            string fusion_id = cand.first;
            std::map<string, vector<locus>> breakpoints;

            bool is_forward = true;
            for(const auto &ff : {cand.second.forward, cand.second.backward}){//, cand.second.no_first, cand.second.multi_first}){
                for(const auto &fus : ff){
                    for(const auto &bp_pair : fus.get_breakpoints(is_forward)){

                        bp_file << fus.read_id <<"\t" << fusion_id << "\t" << bp_pair.first << "\t" << bp_pair.second << "\n";
                        breakpoints[bp_pair.first].push_back(bp_pair.second);
                    }
                }
                is_forward = false;
            }
        }
        bp_file.close();
        return 0;  
    }
       
} 
