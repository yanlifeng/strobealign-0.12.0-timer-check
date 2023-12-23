// ssw_cpp.cpp
// Created by Wan-Ping Lee
// Last revision by Mengyao Zhao on 2023-Apr-21

#include "ssw_cpp.h"
#include "ssw.h"

#include <fstream>  // 包含必要的头文件
#include <cstring>  // 如果您还没有包含此头文件，请添加它
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>
#include <map>
#include <string>
#include <chrono>  // 如果您使用了时间相关的功能

#include <mutex>
#include <condition_variable>


namespace {

    static const int8_t kBaseTranslation[128] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        //   A     C            G
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        //             T
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        //   a     c            g
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        //             t
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };

    void BuildSwScoreMatrix(const uint8_t& match_score,
            const uint8_t& mismatch_penalty,
            int8_t* matrix) {

        // The score matrix looks like
        //                 // A,  C,  G,  T,  N
        //  score_matrix_ = { 2, -2, -2, -2, -2, // A
        //                   -2,  2, -2, -2, -2, // C
        //                   -2, -2,  2, -2, -2, // G
        //                   -2, -2, -2,  2, -2, // T
        //                   -2, -2, -2, -2, -2};// N

        int id = 0;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                matrix[id] = ((i == j) ? match_score : static_cast<int8_t>(-mismatch_penalty));
                ++id;
            }
            matrix[id] = static_cast<int8_t>(-mismatch_penalty); // For N
            ++id;
        }

        for (int i = 0; i < 5; ++i)
            matrix[id++] = static_cast<int8_t>(-mismatch_penalty); // For N

    }

    void ConvertAlignment(const s_align& s_al,
            const int& query_len,
            StripedSmithWaterman::Alignment* al) {
        al->sw_score           = s_al.score1;
        al->sw_score_next_best = s_al.score2;
        al->ref_begin          = s_al.ref_begin1;
        al->ref_end            = s_al.ref_end1;
        al->query_begin        = s_al.read_begin1;
        al->query_end          = s_al.read_end1;
        al->ref_end_next_best  = s_al.ref_end2;

        al->cigar.clear();
        al->cigar_string.clear();

        if (s_al.cigarLen > 0) {
            std::ostringstream cigar_string;
            if (al->query_begin > 0) {
                uint32_t cigar = to_cigar_int(al->query_begin, 'S');
                al->cigar.push_back(cigar);
                cigar_string << al->query_begin << 'S';
            }

            for (int i = 0; i < s_al.cigarLen; ++i) {
                al->cigar.push_back(s_al.cigar[i]);
                cigar_string << cigar_int_to_len(s_al.cigar[i]) << cigar_int_to_op(s_al.cigar[i]);
            }

            int end = query_len - al->query_end - 1;
            if (end > 0) {
                uint32_t cigar = to_cigar_int(end, 'S');
                al->cigar.push_back(cigar);
                cigar_string << end << 'S';
            }

            al->cigar_string = cigar_string.str();
        } // end if
    }

    // @Function:
    //     Calculate the length of the previous cigar operator
    //     and store it in new_cigar and new_cigar_string.
    //     Clean up in_M (false), in_X (false), length_M (0), and length_X(0).
    void CleanPreviousMOperator(
            bool* in_M,
            bool* in_X,
            uint32_t* length_M,
            uint32_t* length_X,
            std::vector<uint32_t>* new_cigar,
            std::ostringstream* new_cigar_string) {
        if (*in_M) {
            //fprintf(stderr, "inM %d\n", *length_M);
            uint32_t match = to_cigar_int(*length_M, '=');
            new_cigar->push_back(match);
            (*new_cigar_string) << *length_M << '=';
        } else if (*in_X){ //in_X
            //fprintf(stderr, "inX %d\n", *length_X);
            uint32_t match = to_cigar_int(*length_X, 'X');
            new_cigar->push_back(match);
            (*new_cigar_string) << *length_X << 'X';
        }

        // Clean up
        *in_M = false;
        *in_X = false;
        *length_M = 0;
        *length_X = 0;
    }

    // @Function:
    //     1. Calculate the number of mismatches.
    //     2. Modify the cigar string (sample):
    // @Return:
    //     The number of mismatches.
    int CalculateNumberMismatchOnly(
            StripedSmithWaterman::Alignment* al,
            int8_t const *ref,
            int8_t const *query,
            const int& query_len) {

        ref   += al->ref_begin;
        query += al->query_begin;
        int mismatch_length = 0;

        //fprintf(stderr, "old str %s\n", al->cigar_string.c_str());

        std::vector<uint32_t> new_cigar;
        std::ostringstream new_cigar_string;

        //if (al->query_begin > 0) {
        //	uint32_t cigar = to_cigar_int(al->query_begin, 'S');
        //	new_cigar.push_back(cigar);
        //	new_cigar_string << al->query_begin << 'S';
        //}


        for (unsigned int i = 0; i < al->cigar.size(); ++i) {
            char op = cigar_int_to_op(al->cigar[i]);
            uint32_t length = cigar_int_to_len(al->cigar[i]);
            if (op == 'M') {
                uint32_t match = to_cigar_int(length, '=');
                new_cigar.push_back(match);
                new_cigar_string << length << '=';
            } else if (op == 'S' || op == 'X' || op == 'I' || op == 'D') {
                uint32_t match = to_cigar_int(length, op);
                new_cigar.push_back(match);
                new_cigar_string << length << op;
                if(op != 'S') mismatch_length += length;
            } else {
                fprintf(stderr, "GGGGGGGGGGGGGG");
                exit(0);
            }
        }


        //int end = query_len - al->query_end - 1;
        //if (end > 0) {
        //	uint32_t cigar = to_cigar_int(end, 'S');
        //	new_cigar.push_back(cigar);
        //	new_cigar_string << end << 'S';
        //}

        al->cigar_string.clear();
        al->cigar.clear();
        al->cigar_string = new_cigar_string.str();
        al->cigar = new_cigar;
        //fprintf(stderr, "new str %s\n", al->cigar_string.c_str());

        return mismatch_length;
    }


    // @Function:
    //     1. Calculate the number of mismatches.
    //     2. Modify the cigar string:
    //         differentiate matches (M) and mismatches(X).
    // @Return:
    //     The number of mismatches.
    int CalculateNumberMismatch(
            StripedSmithWaterman::Alignment* al,
            int8_t const *ref,
            int8_t const *query,
            const int& query_len) {

        ref   += al->ref_begin;
        query += al->query_begin;
        int mismatch_length = 0;

        //fprintf(stderr, "old str %s\n", al->cigar_string.c_str());

        std::vector<uint32_t> new_cigar;
        std::ostringstream new_cigar_string;

        if (al->query_begin > 0) {
            uint32_t cigar = to_cigar_int(al->query_begin, 'S');
            new_cigar.push_back(cigar);
            new_cigar_string << al->query_begin << 'S';
        }

        bool in_M = false; // the previous is match
        bool in_X = false; // the previous is mismatch
        uint32_t length_M = 0;
        uint32_t length_X = 0;

        for (unsigned int i = 0; i < al->cigar.size(); ++i) {
            char op = cigar_int_to_op(al->cigar[i]);
            uint32_t length = cigar_int_to_len(al->cigar[i]);
            //fprintf(stderr, "%c %d\n", op, length);
            if (op == 'M') {
                for (uint32_t j = 0; j < length; ++j) {
                    //fprintf(stderr, "%d %d %d\n", j, *ref, *query);

                    if (*ref != *query) {
                        ++mismatch_length;
                        if (in_M) { // the previous is match; however the current one is mismatche
                            uint32_t match = to_cigar_int(length_M, '=');
                            new_cigar.push_back(match);
                            new_cigar_string << length_M << '=';
                            //fprintf(stderr, "inM %d\n", length_M);
                        }
                        length_M = 0;
                        ++length_X;
                        in_M = false;
                        in_X = true;
                    } else { // *ref == *query
                        if (in_X) { // the previous is mismatch; however the current one is matche
                            uint32_t match = to_cigar_int(length_X, 'X');
                            new_cigar.push_back(match);
                            new_cigar_string << length_X << 'X';
                            //fprintf(stderr, "inX %d\n", length_X);
                        }
                        ++length_M;
                        length_X = 0;
                        in_M = true;
                        in_X = false;
                    } // end of if (*ref != *query)
                    ++ref;
                    ++query;
                }
            } else if (op == 'I') {
                query += length;
                mismatch_length += length;
                CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);
                new_cigar.push_back(al->cigar[i]);
                new_cigar_string << length << 'I';
            } else if (op == 'D') {
                ref += length;
                mismatch_length += length;
                CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);
                new_cigar.push_back(al->cigar[i]);
                new_cigar_string << length << 'D';
            }
        }

        CleanPreviousMOperator(&in_M, &in_X, &length_M, &length_X, &new_cigar, &new_cigar_string);

        int end = query_len - al->query_end - 1;
        if (end > 0) {
            uint32_t cigar = to_cigar_int(end, 'S');
            new_cigar.push_back(cigar);
            new_cigar_string << end << 'S';
        }

        al->cigar_string.clear();
        al->cigar.clear();
        al->cigar_string = new_cigar_string.str();
        al->cigar = new_cigar;
        //fprintf(stderr, "new str %s\n", al->cigar_string.c_str());

        return mismatch_length;
    }

    void SetFlag(const StripedSmithWaterman::Filter& filter, uint8_t* flag) {
        if (filter.report_begin_position) *flag |= 0x08;
        if (filter.report_cigar) *flag |= 0x0f;
    }

    // http://www.cplusplus.com/faq/sequences/arrays/sizeof-array/#cpp
    template <typename T, size_t N>
        inline size_t SizeOfArray( const T(&)[ N ] )
        {
            return N;
        }

} // namespace



namespace StripedSmithWaterman {

    Aligner::Aligner(void)
        : score_matrix_(NULL)
          , score_matrix_size_(5)
          , translation_matrix_(NULL)
          , match_score_(2)
          , mismatch_penalty_(2)
          , gap_opening_penalty_(3)
          , gap_extending_penalty_(1)
          , translated_reference_(NULL)
          , reference_length_(0)
    {
        BuildDefaultMatrix();
    }

    Aligner::Aligner(
            const uint8_t& match_score,
            const uint8_t& mismatch_penalty,
            const uint8_t& gap_opening_penalty,
            const uint8_t& gap_extending_penalty)

        : score_matrix_(NULL)
        , score_matrix_size_(5)
        , translation_matrix_(NULL)
        , match_score_(match_score)
        , mismatch_penalty_(mismatch_penalty)
        , gap_opening_penalty_(gap_opening_penalty)
        , gap_extending_penalty_(gap_extending_penalty)
        , translated_reference_(NULL)
        , reference_length_(0)
        {
            BuildDefaultMatrix();
        }

    Aligner::Aligner(const int8_t* score_matrix,
            const int&    score_matrix_size,
            const int8_t* translation_matrix,
            const int&    translation_matrix_size)

        : score_matrix_(NULL)
        , score_matrix_size_(score_matrix_size)
        , translation_matrix_(NULL)
        , match_score_(2)
        , mismatch_penalty_(2)
        , gap_opening_penalty_(3)
        , gap_extending_penalty_(1)
        , translated_reference_(NULL)
        , reference_length_(0)
        {
            score_matrix_ = new int8_t[score_matrix_size_ * score_matrix_size_];
            memcpy(score_matrix_, score_matrix, sizeof(int8_t) * score_matrix_size_ * score_matrix_size_);
            translation_matrix_ = new int8_t[translation_matrix_size];
            memcpy(translation_matrix_, translation_matrix, sizeof(int8_t) * translation_matrix_size);
        }


    Aligner::~Aligner(void){
        Clear();
    }

    int Aligner::SetReferenceSequence(const char* seq, const int& length) {

        int len = 0;
        if (translation_matrix_) {
            // delete the current buffer
            CleanReferenceSequence();
            // allocate a new buffer
            translated_reference_ = new int8_t[length];

            len = TranslateBase(seq, length, translated_reference_);
        } else {
            // nothing
        }

        reference_length_ = len;
        return len;
    }

    int Aligner::TranslateBase(const char* bases, const int& length,
            int8_t* translated) const {

        const char* ptr = bases;
        int len = 0;
        for (int i = 0; i < length; ++i) {
            translated[i] = translation_matrix_[(int) *ptr];
            ++ptr;
            ++len;
        }

        return len;
    }


    uint16_t Aligner::Align(const char* query, const Filter& filter,
            Alignment* alignment, const int32_t maskLen) const
    {
        if (!translation_matrix_) return false;
        if (reference_length_ == 0) return false;

        int query_len = strlen(query);
        if (query_len == 0) return false;
        int8_t* translated_query = new int8_t[query_len];
        TranslateBase(query, query_len, translated_query);

        const int8_t score_size = 2;
        s_profile* profile = ssw_init(translated_query, query_len, score_matrix_,
                score_matrix_size_, score_size);

        uint8_t flag = 0;
        SetFlag(filter, &flag);
        s_align* s_al = ssw_align(profile, translated_reference_, reference_length_,
                static_cast<int>(gap_opening_penalty_),
                static_cast<int>(gap_extending_penalty_),
                flag, filter.score_filter, filter.distance_filter, maskLen);

        alignment->Clear();
        ConvertAlignment(*s_al, query_len, alignment);
        alignment->mismatches = CalculateNumberMismatch(&*alignment, translated_reference_, translated_query, query_len);
        uint16_t align_flag = s_al->flag;

        // Free memory
        delete [] translated_query;
        align_destroy(s_al);
        init_destroy(profile);

        return align_flag;
    }


#include <sys/time.h>


    inline double GetTime() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
    }

    int read_last_two_lines(const std::string& filename, std::string &name) {
        std::ifstream file(filename, std::ios::ate); // Open file for reading at the end
        if (!file.is_open()) {
            std::cerr << "Unable to open file: " << filename << std::endl;
            return 0;
        }

        std::string lastLine, secondLastLine;
        long long fileSize = file.tellg();

        if (fileSize <= 0) {
            std::cerr << "File is empty or an error occurred." << std::endl;
            return 0;
        }

        char ch;
        int newLinesCount = 0;

        // Iterate backwards through the file
        for (long long i = fileSize - 1; i >= 0; --i) {
            file.seekg(i);
            file.get(ch);

            if (ch == '\n') {
                newLinesCount++;
                if (newLinesCount > 2) {
                    break;
                }
            }
        }

        // Reading the lines after the second last '\n' found
        if (newLinesCount > 2) {
            getline(file, secondLastLine);
            getline(file, lastLine);
        } else if (newLinesCount == 2) {
            file.seekg(0);
            getline(file, lastLine);
            getline(file, lastLine);
        } else if(newLinesCount == 1) {
            file.seekg(0);
            getline(file, lastLine);
        } else {
            std::cerr << "Error count " << newLinesCount << std::endl;
            return 0;
        }

        //printf("secondLastLine : %s\n", secondLastLine.c_str());
        //printf("lastLine : %s\n", lastLine.c_str());


        file.close();

        if (lastLine.substr(0, 5) == "name:") {
            name = lastLine.substr(6);
            return 0;
        } if (lastLine.substr(0, 4) == "ame:") {
            name = lastLine.substr(5);
            return 0;
        } else {
            std::size_t separatorPos = secondLastLine.find('_');
            if (separatorPos != std::string::npos) {
                name = secondLastLine.substr(1, separatorPos - 1);
                int id = std::stoi(secondLastLine.substr(separatorPos + 1));
                return id + 1;
            } else {
                printf("secondLastLine : %s\n", secondLastLine.c_str());
                printf("lastLine : %s\n", lastLine.c_str());
                std::cerr << "read name error" << std::endl; 
            }
        }
    }

    struct Data {
        int value1;
        int value2;
        int value3;
        int value4;
        int value5;
        std::string value6;

    };
    std::mutex dataMapMutex;
    std::condition_variable dataMapCreated;
    bool isDataMapReady = false;
    std::map<std::pair<int, int>, Data> dataMap; // 全局变量


    int calculate_cigar_length(const char *cigar) {
        int length = 0;
        int current_length = 0;

        while (*cigar != '\0') {
            if (isdigit(*cigar)) {
                current_length = current_length * 10 + (*cigar - '0');
            } else {
                if (*cigar == 'M' || *cigar == 'I' || *cigar == 'X' || *cigar == 'S') {
                    length += current_length;
                }
                current_length = 0;
            }
            cigar++;
        }

        return length;
    }
    uint16_t Aligner::Align(const char* query, const char* ref, const int& ref_len,
            const Filter& filter, Alignment* alignment, const int32_t maskLen) const
    {
        static thread_local double init_time = 0;
        static thread_local double ssw_time = 0;
        static thread_local double convert_time = 0;
        static thread_local double cal_time = 0;
        static thread_local double tot_time = 0;
        static thread_local int cnt = 0;
        double t0;


        // 打开文件，使用append模式
        std::thread::id mainThreadId = std::this_thread::get_id();
        std::ostringstream ss;
        ss << mainThreadId;
        std::string threadIdStr = ss.str();
        unsigned long mainThreadIdValue = std::stoul(ss.str());
        int has_v6 = 1;


        //this ccc
#ifdef RRUUNN
        if(cnt == 0) {
#else
        if(cnt == -1) {
#endif

            std::unique_lock<std::mutex> lock(dataMapMutex);
            if (!isDataMapReady) {
                fprintf(stderr, "%ld begin\n", mainThreadIdValue);
                // 0号线程创建 dataMap
                // ...创建 dataMap 的代码...
                std::ifstream file("merge.res"); // 替换为您的文件名
                std::string line;

                while (getline(file, line)) {
                    std::istringstream iss(line);
                    std::string prefix1, prefix2, srr;
                    int v1, v2, v3, v4, v5;
                    std::string v6 = "";

                    if(has_v6) {
                        if (!(iss >> prefix1 >> prefix2 >> srr >> v1 >> v2 >> v3 >> v4 >> v5 >> v6)) {
                            continue; // 解析错误或行格式不匹配
                        }
                    } else {
                        if (!(iss >> prefix1 >> prefix2 >> srr >> v1 >> v2 >> v3 >> v4 >> v5)) {
                            continue; // 解析错误或行格式不匹配
                        }
                    }


                    // 提取SRR或RR后的数字
                    size_t pos = srr.find('.');
                    size_t underscore = srr.find('_');
                    if (pos == std::string::npos || underscore == std::string::npos || pos >= underscore) {
                        std::cerr << "Format error " << srr << std::endl;
                        continue; // 格式错误
                    }
                    int num1 = std::stoi(srr.substr(pos + 1, underscore - pos - 1));
                    int num2 = std::stoi(srr.substr(underscore + 1));
                    //std::cerr << srr << " " << num1 << " " << num2 << std::endl;

                    if(has_v6) {
                        std::string cigar_substring = v6.substr(6);
                        const char *cigar_str = cigar_substring.c_str();
                        int seq_length = calculate_cigar_length(cigar_str);
                        if(seq_length != v4 - v2 + 1) {
                            v1 = 0; v2 = 0; v3 = 0; v4 = 0; v5 = 0;
                            //fprintf(stderr, "oh GG len %d %d\n", seq_length, v4 - v2 + 1);
                        }
                    }
                    if(v4 < v2 || v5 < v3 || v2 < 0 || v3 < 0 || v4 < 0 || v5 < 0 || v1 == 0) {
                        //fprintf(stderr, "oh GG v4 v2 v5 v3 value gg\n");
                        v1 = 0; v2 = 0; v3 = 0; v4 = 0; v5 = 0;
                    }
                    // 存储数据
                    dataMap[{num1, num2}] = {v1, v2, v3, v4, v5, v6};
                }

                isDataMapReady = true;
                dataMapCreated.notify_all(); // 通知其他线程
                fprintf(stderr, "%ld done\n", mainThreadIdValue);
            } else {
                // 其他线程等待 dataMap 创建完成
                dataMapCreated.wait(lock, []{ return isDataMapReady; });
            }
            lock.unlock();

        }
        std::string queryFilename = "query_sequences_" + threadIdStr + ".fasta";
        std::string refFilename = "ref_sequences_" + threadIdStr + ".fasta";
        //fprintf(stderr, "queryFilename %s\n", queryFilename.c_str());

        int cntt = 0;

        //fprintf(stderr, "=== [debug ssw] %d\n", cntt++);
        std::string read_name;
        int read_id = -1;

        std::ofstream query_file(queryFilename, std::ios::app);
        std::ofstream ref_file(refFilename, std::ios::app);
        read_id = read_last_two_lines(queryFilename, read_name);
        //fprintf(stderr, "read name: %s\n", read_name.c_str());
        //fprintf(stderr, "read_id %d\n", read_id);


        if (query_file.is_open() && ref_file.is_open()) {
            // 写入query序列
            query_file << ">" << read_name << "_" << read_id << "\n";  // 写入序列标识符
            query_file << query << "\n";  // 写入序列数据

            // 写入ref序列
            ref_file << ">" << read_name << "_" << read_id << "\n";       // 写入序列标识符
            ref_file << ref << "\n";  // 写入序列数据
        } else {
            // 如果文件无法打开，输出错误信息
            std::cerr << "Error: Unable to open the FASTA files for writing." << std::endl;
        }

        // 关闭文件
        query_file.close();
        ref_file.close();

        //fprintf(stderr, "=== [debug ssw] %d\n", cntt++);
        int pp_pos = read_name.find('.');
        int name_id = std::stoi(read_name.substr(pp_pos + 1));
        auto key = std::make_pair(name_id, read_id);
        Data data = {0, 0, 0, 0, 0, ""};
        //this ccc
#ifdef RRUUNN
        if (dataMap.find(key) == dataMap.end()) {
            std::cerr << "GG no " << name_id << " " << read_id << std::endl;
        } else {
            data = dataMap[key];
        }
#endif
        //std::cerr << "Values: " << data.value1 << " " << data.value2 << " " << data.value3 << " "
        //    << data.value4 << " " << data.value5 << std::endl;



        //fprintf(stderr, "=== [debug ssw] %d\n", cntt++);

        double tt0 = GetTime();
        if (!translation_matrix_) return false;

        int query_len = strlen(query);
        if (query_len == 0) return false;
        int8_t* translated_query = new int8_t[query_len];
        TranslateBase(query, query_len, translated_query);

        // calculate the valid length
        int valid_ref_len = ref_len;
        int8_t* translated_ref = new int8_t[valid_ref_len];
        TranslateBase(ref, valid_ref_len, translated_ref);



        const int8_t score_size = 2;
        t0 = GetTime();
        s_profile* profile = ssw_init(translated_query, query_len, score_matrix_,
                score_matrix_size_, score_size);
        init_time += GetTime() - t0;

        uint8_t flag = 0;
        SetFlag(filter, &flag);
        t0 = GetTime();

        s_align* s_al; 
        //this ccc
#ifdef RRUUNN
        if(data.value1 + data.value2 + data.value3 + data.value4 + data.value5 == 0) {
#else
        if(1) {
#endif
            s_al = ssw_align(profile, translated_ref, valid_ref_len,
                    static_cast<int>(gap_opening_penalty_),
                    static_cast<int>(gap_extending_penalty_),
                    flag, filter.score_filter, filter.distance_filter, maskLen);
        } else {
            s_al = my_ssw_align(profile, translated_ref, valid_ref_len,
                    static_cast<int>(gap_opening_penalty_),
                    static_cast<int>(gap_extending_penalty_),
                    flag, filter.score_filter, filter.distance_filter, maskLen, data.value1, data.value2, data.value3, data.value4, data.value5, data.value6.c_str());
        }
        //fprintf(stderr, " [info %d] %d %d %d %d %d %d %d %d\n", cnt, query_len, ref_len, static_cast<int>(gap_opening_penalty_), static_cast<int>(gap_extending_penalty_), flag, filter.score_filter, filter.distance_filter, maskLen);
        //fprintf(stderr, " [res %d] score:%d, query_begin:%d, ref_begin:%d, query_end:%d, ref_end:%d\n", cnt, s_al->score1, s_al->read_begin1, s_al->ref_begin1, s_al->read_end1, s_al->ref_end1);
        ssw_time += GetTime() - t0;

        alignment->Clear();
        t0 = GetTime();
        ConvertAlignment(*s_al, query_len, alignment);
        //fprintf(stderr, "s2 : %s\n", alignment->cigar_string.c_str());
        convert_time += GetTime() - t0;
        t0 = GetTime();
        if(data.value6.length()) {
            alignment->mismatches = CalculateNumberMismatchOnly(&*alignment, translated_ref, translated_query, query_len);
        } else {
            alignment->mismatches = CalculateNumberMismatch(&*alignment, translated_ref, translated_query, query_len);
        }
        //fprintf(stderr, "s3 : %s\n", alignment->cigar_string.c_str());
        cal_time += GetTime() - t0;
        uint16_t align_flag = s_al->flag;

        // Free memory
        delete [] translated_query;
        delete [] translated_ref;
        align_destroy(s_al);
        init_destroy(profile);
        tot_time += GetTime() - tt0;
        //fprintf(stderr, "=== [debug ssw] %d\n", cntt++);
        //if(cnt == 96203) {
        //    exit(0);
        //}
        cnt++;

        if(cnt % 10000 == 0) {
            std::thread::id mainThreadId = std::this_thread::get_id();
            std::ostringstream ss;
            ss << mainThreadId;
            unsigned long mainThreadIdValue = std::stoul(ss.str());
            fprintf(stderr, "==== [%lld] tot:%lf init:%lf ssw:%lf convert:%lf cal:%lf\n", mainThreadIdValue, tot_time, init_time, ssw_time, convert_time, cal_time);
        }


        //fprintf(stderr, "=== [debug ssw] %d\n", cntt++);
        return align_flag;
        }

        void Aligner::Clear(void) {
            ClearMatrices();
            CleanReferenceSequence();
        }

        void Aligner::SetAllDefault(void) {
            score_matrix_size_     = 5;
            match_score_           = 2;
            mismatch_penalty_      = 2;
            gap_opening_penalty_   = 3;
            gap_extending_penalty_ = 1;
            reference_length_      = 0;
        }

        bool Aligner::ReBuild(void) {
            if (translation_matrix_) return false;

            SetAllDefault();
            BuildDefaultMatrix();

            return true;
        }

        bool Aligner::ReBuild(
                const uint8_t& match_score,
                const uint8_t& mismatch_penalty,
                const uint8_t& gap_opening_penalty,
                const uint8_t& gap_extending_penalty) {
            if (translation_matrix_) return false;

            SetAllDefault();

            match_score_           = match_score;
            mismatch_penalty_      = mismatch_penalty;
            gap_opening_penalty_   = gap_opening_penalty;
            gap_extending_penalty_ = gap_extending_penalty;

            BuildDefaultMatrix();

            return true;
        }

        bool Aligner::ReBuild(
                const int8_t* score_matrix,
                const int&    score_matrix_size,
                const int8_t* translation_matrix,
                const int&    translation_matrix_size) {

            ClearMatrices();
            score_matrix_ = new int8_t[score_matrix_size_ * score_matrix_size_];
            memcpy(score_matrix_, score_matrix, sizeof(int8_t) * score_matrix_size_ * score_matrix_size_);
            translation_matrix_ = new int8_t[translation_matrix_size];
            memcpy(translation_matrix_, translation_matrix, sizeof(int8_t) * translation_matrix_size);

            return true;
        }

        void Aligner::BuildDefaultMatrix(void) {
            ClearMatrices();
            score_matrix_ = new int8_t[score_matrix_size_ * score_matrix_size_];
            BuildSwScoreMatrix(match_score_, mismatch_penalty_, score_matrix_);
            translation_matrix_ = new int8_t[SizeOfArray(kBaseTranslation)];
            memcpy(translation_matrix_, kBaseTranslation, sizeof(int8_t) * SizeOfArray(kBaseTranslation));
        }

        void Aligner::ClearMatrices(void) {
            delete [] score_matrix_;
            score_matrix_ = NULL;

            delete [] translation_matrix_;
            translation_matrix_ = NULL;
        }
        } // namespace StripedSmithWaterman
