#include <iostream>
#include <list>

int main() {
    std::list<int> myList = {1, 2, 3, 4, 5};

    for (auto it = myList.begin(); it != myList.end(); /* no increment here */) {
        if (*it % 2 == 0) { // Condition to delete an element
            it = myList.erase(it); // Erase element and update iterator
        } else {
            ++it; // Only increment iterator if no deletion
        }

        // Example condition to add an element
        if (*it == 3) {
            myList.insert(it, 10); // Insert before the current element
        }
    }

    // Printing the list
    for (int n : myList) {
        std::cout << n << ' ';
    }

    return 0;
}