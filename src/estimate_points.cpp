// #include <Eigen/Dense>

#include <iostream>


#include <string>
#include <vector>
#include <map>


#include <chrono>


#include "file_reading.cpp"
#include "est_algo.cpp"




int main(int argc, char **argv)
{

    // std::string txt_path = "../Task-Points/3-lines.txt";

    std::string txt_path = argv[1];
    std::vector<std::vector<double>> lines;
    int n_points;
    read_txt(txt_path, lines, n_points);

    auto estimation_start = std::chrono::high_resolution_clock::now();
    estimate_points_greedy(lines,n_points);
    auto estimation_end = std::chrono::high_resolution_clock::now();

    auto estimation_duration = std::chrono::duration_cast<std::chrono::milliseconds>(estimation_end - estimation_start);

    std::cout << "estimation_duration: " << estimation_duration.count() << " milliseconds" << std::endl;


    return 0;
}