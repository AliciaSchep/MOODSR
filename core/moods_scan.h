#ifndef MOODS_SCAN_H
#define MOODS_SCAN_H

#include <string>

#include "moods.h"

namespace MOODS { namespace scan{
    
    struct match
    {
        size_t pos;
        double score;
    };
    
    std::vector< std::vector< match> > scan_dna(const std::string& seq,
                                                const std::vector<score_matrix>& matrices,
                                                const std::vector<double>& bg,
                                                const std::vector<double> thresholds,
                                                unsigned int window_size );

    std::vector< std::vector< match> > scan(const std::string& seq,
                                            const std::vector<score_matrix>& matrices,
                                            const std::vector<double>& bg,
                                            const std::vector<double> thresholds,
                                            unsigned int window_size,
                                            const std::vector<std::string>& alphabet);
    // unsigned int scan_dna(const std::string& seq, const std::vector<score_matrix>& matrices, const std::vector<double>& bg, const std::vector<double> thresholds, unsigned int window_size );


    
}}



#endif