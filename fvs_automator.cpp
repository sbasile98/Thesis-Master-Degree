#include <bits/stdc++.h>
#include <dirent.h>

using namespace std;

const char *get_filename_ext (const char *filename) {
    const char *dot = strrchr(filename, '.');
    if (!dot || dot == filename) return "";
    return dot + 1;
}

int main(int argc, char *argv[]) {
    DIR *dir;
    struct dirent *ent;
	
	double delta = 1;
	double epsilon = 0.5;
	double lambda = 0.2;
	double starting_mutation = 30.0;
    
    string directory(argv[1]);
    if (directory.back() != '/') directory += "/";
    
    if ((dir = opendir(argv[1])) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (strcmp(get_filename_ext(ent->d_name), "fvs") == 0) {
                string filename(ent->d_name);
				system(("./fvs_ml.out -f " + directory + filename + " -e " + to_string(epsilon) + " -d " + to_string(delta) + " -l " + to_string(lambda) + " -m " + to_string(starting_mutation)).c_str());
            }
        }
        closedir(dir);
    }
    else {
        perror("");
        return EXIT_FAILURE;
    }
    
    return 0;
}
