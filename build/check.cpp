#include <fstream>
#include <iostream>
#include <set>
#include <string>

bool checkFile(const std::string& filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "无法打开文件：" << filename << std::endl;
		return false;
	}

	std::string line;
	std::set<std::string> queries;
	bool inNameSection = false;

	while (getline(file, line)) {
		if (line.substr(0, 5) == "name:") {
			if (queries.size() > 4) {
                printf("GG %d\n", queries.size());
                //printf("GG %s\n", line.c_str());
                //for(auto S : queries) printf("%s\n", S.c_str());
				return false; // 发现超过两种查询序列
			}
			queries.clear();
		} else if (line[0] == '>') {
		} else{
			queries.insert(line);
		}
	}

	return queries.size() <= 2;
}

int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cerr << "用法：" << argv[0] << " <filename>" << std::endl;
		return 1;
	}

	std::string filename = argv[1];
	if (checkFile(filename)) {
		std::cout << "文件格式正确" << std::endl;
	} else {
		std::cout << "文件格式不正确" << std::endl;
	}

	return 0;
}

