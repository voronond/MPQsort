#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include <aqsort.h>

int main()
{
    std::vector<uint32_t> rows {   0,   3,   1,   2,   0,   3 };
    std::vector<uint32_t> cols {   0,   0,   1,   2,   3,   3 };
    std::vector<double>   vals { 1.0, 5.0, 3.0, 4.0, 2.0, 6.0 };

    for (int i = 0; i < rows.size(); i++) 
        std::cout << "(" << rows[i] << ", " << cols[i] << ") = " << vals[i] << std::endl;
    std::cout << std::endl;

    auto comp = [&rows, &cols] (std::size_t i, std::size_t j) /* -> bool */ {
        if (rows[i] < rows[j])
            return true;
        if ((rows[i] == rows[j]) && (cols[i] < cols[j]))
            return true;
        return false;
    };

    auto swap = [&rows, &cols, &vals] (std::size_t i, std::size_t j) {
        std::swap(rows[i], rows[j]);
        std::swap(cols[i], cols[j]);
        std::swap(vals[i], vals[j]);
    };

    aqsort::sort(rows.size(), &comp, &swap);

    for (int i = 0; i < rows.size(); i++) 
        std::cout << "(" << rows[i] << ", " << cols[i] << ") = " << vals[i] << std::endl;
}
