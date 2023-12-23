#include <iostream>
#include <fstream>
#include <string>
#include <vector>

void read_last_two_lines(const std::string& filename) {
    std::ifstream file(filename, std::ios::ate); // Open file for reading at the end
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return;
    }

    std::string lastLine, secondLastLine;
    long long fileSize = file.tellg();

    if (fileSize <= 0) {
        std::cerr << "File is empty or an error occurred." << std::endl;
        return;
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
    printf("newLinesCount %d\n", newLinesCount);

    // Reading the lines after the second last '\n' found
    if (newLinesCount > 2) {
        getline(file, secondLastLine);
        getline(file, lastLine);
    } else if (newLinesCount <= 2) {
        //file.seekg(0);
        getline(file, lastLine);
    } else {
        std::cerr << "Error count " << newLinesCount << std::endl;
        return;
    }
    //if(newLinesCount == 1) std::swap(secondLastLine, lastLine);

    std::cout << "Second Last Line: " << secondLastLine << std::endl;
    std::cout << "Last Line: " << lastLine << std::endl;

    file.close();
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    read_last_two_lines(filename);
    return 0;
}

