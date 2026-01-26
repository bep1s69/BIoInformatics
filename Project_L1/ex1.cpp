#include <iostream>
#include <set>
#include <string>

int main() {
    std::string seq = "A T T G C C C C G A A T";
    std::set<char> alphabet;
    
        for (char c : seq) {
        if (c != ' ') {
            alphabet.insert(c);
        }
    }
    
        std::cout << "Alphabet of the sequence: ";
    for (char c : alphabet) {
        std::cout << c << " ";
    }
    std::cout << "\n";
    
    return 0;
}