

#include <fstream>
#include <cassert>

/*
    splits string by char 
*/
std::vector<std::string> split_by_char(std::string input_line, char splitter)
{
    std::vector<std::string> splitted_input;

    std::string input_cur;
    input_cur.clear();
    for(auto x : input_line)
    {
        // std::cout << "x: " << x << "input_cur: " << input_cur << std::endl;
        if(x == splitter)
        {
            splitted_input.push_back(input_cur);
            input_cur.clear();
        }
        else
        {
            input_cur += x;
        }
    }

    splitted_input.push_back(input_cur);

    return splitted_input;
}

/*
    reads txt from input_name path line by line
    and writes everything into lines vector of vectors, where
    outer idx - line index
    inner idx - {x_1,y_1,x_2,y_2} , represent 2 points building a line
*/
void read_txt(const std::string &input_name, std::vector<std::vector<double>> &lines, int &n_points)
{
    std::ifstream in_file;
    std::string in_line;


    // in_file.open("../Task-Points/3-lines.txt");
    in_file.open(input_name);

    bool first_line = true;
    int n_lines;

    std::vector<std::string> input_splitted;
    // std::vector<std::vector<double>> lines;

    // read txt line by line
    while(std::getline(in_file,in_line))
    {
        input_splitted = split_by_char(in_line,' ');
        if(first_line)
        {
            first_line = false;
            n_lines = std::stod(input_splitted[0]);
            n_points = std::stod(input_splitted[1]);
        }
        else
        {
            double x_1 = std::stod(input_splitted[0]);
            double y_1 = std::stod(input_splitted[1]);
            double x_2 = std::stod(input_splitted[2]);
            double y_2 = std::stod(input_splitted[3]);

            std::vector<double> line = {x_1,y_1,x_2,y_2};

            lines.push_back(line);
        }
    }

    assert(n_lines == lines.size());

    std::cout << "n_lines: " << n_lines << std::endl;
    std::cout << "n_points: " << n_points << std::endl;

}