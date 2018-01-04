/**
 * Inverse word search program written by Nicholas Roper in 2017
 * for Data Structures HW6 contest at RPI, for which I was named
 * among the "Most Fastest" and "Most Correct". For the 291 submitted
 * programs, of those that correctly solved every test case, this
 * program was one of 2 whose runtimes were rounded down to 0.0s
 * for each test case except the last, and it had the lowest
 * overall time taken to run each test case.
 * See http://www.cs.rpi.edu/academics/courses/spring17/ds/hw6_contest_results.php
 * for more precise details.
 * Finds one/all grids of a certain width and height that contain
 * every word in an input list of included words and contain no
 * words from another input list of excluded words.
 *
 * Input file format is:
 
width height
+ include_word1
+ include_word2
+ include_word3
...
- exclude_word1
- exclude_word2
- exclude_word3
...

 * where width and height are ints greater than 0, include_word refers to a word that must be
 * somewhere in the puzzle, and exclude_word refers to a word that must not be in the puzzle.
 * All include_words and exclude_words must consist of only chars a-z;
 */


#include <iostream>
// std::cerr
#include <fstream>
// File I/O
#include <unordered_set>
// unordered_set<std::string> is used to hold result grids.
#include <algorithm>
// std::reverse, std::iter_swap

// Basic Node struct for linked list.
struct Node {
	std::string value;
	Node* ptr;
};

// Deallocates linked list. Note that node must be the head of the linked list.
void del_linked_list(Node* node) {
	Node* temp_node;
	while(node) {
		temp_node = node;
		node = node->ptr;
		delete temp_node;
	}
}

/**
 * Used for the first set of checks for an excluded word after all included words have been added.
 * Goes through the entire grid, and checks every possible location/rotation that a word could be.
 * Only checks East, South, Southeast, and Southwest, but this is fine because it checks both
 * forwards and reverse simultaneously in each location, so it checks every possible direction.
 * returns true if the word is in the grid, false if it is not.
 */
bool is_word_in_check1(unsigned int width, unsigned int height, const std::string& Grid_str, const std::string& word) {
	unsigned int word_size = word.size();
	if(height == 0) return false;
	if(word_size > height && word_size > width) return false;
	bool is_for = false;
	bool is_rev = false;
	bool is_dir_for;
	bool is_dir_rev;
	for(unsigned int i = 0; i < height; i++) {
		for(unsigned int j = 0; j < width; j++) {
			if((height-i >= word_size || width-j >= word_size) && ((is_for = (Grid_str[i*width+j] == word[0])) || (is_rev = (Grid_str[i*width+j] == word[word_size-1])))) {
				if(word_size == 1) return true;
				if(is_for) is_rev = (Grid_str[i*width+j] == word[word_size-1]);
				//E
				if(j < width-1 && ((is_dir_for = (is_for && Grid_str[i*width+j+1] == word[1])) || ( is_dir_rev = (is_rev && Grid_str[i*width+j+1] == word[word_size-2]))) && width-j >= word_size) {
					if(is_dir_for) is_dir_rev = is_rev && (Grid_str[i*width+j+1] == word[word_size-2]);
					for(unsigned int k = 2; k < word_size; ++k) {
						if(Grid_str[i*width+j+k] != word[k])
							is_dir_for = false;
						if(Grid_str[i*width+j+k] != word[word_size-k-1])
							is_dir_rev = false;
					}
					if(is_dir_for || is_dir_rev) return true;
				}
				//S
				if(i < height-1 && ((is_dir_for = (is_for && Grid_str[(i+1)*width+j] == word[1])) || ( is_dir_rev = (is_rev && Grid_str[(i+1)*width+j] == word[word_size-2]))) && height-i >= word_size) {
					if(is_dir_for) is_dir_rev = is_rev && (Grid_str[(i+1)*width+j] == word[word_size-2]);
					for(unsigned int k = 2; k < word_size; ++k) {
						if(Grid_str[(i+k)*width+j] != word[k])
							is_dir_for = false;
						if(Grid_str[(i+k)*width+j] != word[word_size-k-1])
							is_dir_rev = false;
					}
					if(is_dir_for || is_dir_rev) return true;
				}
				//SE
				if(j < width-1 && i < height-1 && ((is_dir_for = (is_for && Grid_str[(i+1)*width+j+1] == word[1])) || ( is_dir_rev = (is_rev && Grid_str[(i+1)*width+j+1] == word[word_size-2]))) && width-j >= word_size && height-i >= word_size) {
					if(is_dir_for) is_dir_rev = is_rev && (Grid_str[(i+1)*width+j+1] == word[word_size-2]);
					for(unsigned int k = 2; k < word_size; ++k) {
						if(Grid_str[(i+k)*width+j+k] != word[k])
							is_dir_for = false;
						if(Grid_str[(i+k)*width+j+k] != word[word_size-k-1])
							is_dir_rev = false;
					}
					if(is_dir_for || is_dir_rev) return true;
				}
				//SW
				if(j > 0 && i < height-1 && ((is_dir_for = (is_for && Grid_str[(i+1)*width+j-1] == word[1])) || ( is_dir_rev = (is_rev && Grid_str[(i+1)*width+j-1] == word[word_size-2]))) && j >= word_size-1 && height-i >= word_size) {
					if(is_dir_for) is_dir_rev = is_rev && (Grid_str[(i+1)*width+j-1] == word[word_size-2]);
					for(unsigned int k = 2; k < word_size; ++k) {
						if(Grid_str[(i+k)*width+j-k] != word[k])
							is_dir_for = false;
						if(Grid_str[(i+k)*width+j-k] != word[word_size-k-1])
							is_dir_rev = false;
					}
					if(is_dir_for || is_dir_rev) return true;
				}
			}
		}
	}
	return false;
}

/**
 * Used after adding a single letter when filling the grid. Because only one letter is being added,
 * only the positions around that letter need to be checked. let_loc is the location of the letter
 * in the grid, and every direction and postion that includes that location is checked.
 * is_word_in_check1 would give the same result, but using this function takes less time, as less
 * positions in the grid need to be checked.
 */
bool is_word_in_check2(unsigned int width, unsigned int height, const std::string& Grid_str, const std::string& word, unsigned int let_loc) {
	unsigned int word_size = word.size();
	if(height == 0) return false;
	if(word_size > height && word_size > width) return false;
	unsigned int loc_i = let_loc / width;
	unsigned int loc_j = let_loc % width;
	bool is_for;
	bool is_rev;
	unsigned int ind;
	for(unsigned int i = 0; i < word_size; i++) {
		//E
		ind = loc_i*width+loc_j+1-word_size+i;
		if(loc_j+1+i >= word_size && loc_j+i < width && ((is_for = (Grid_str[ind] == word[0])) || (is_rev = (Grid_str[ind] == word[word_size-1])))) {
			if(is_for) is_rev = Grid_str[ind] == word[word_size-1];
			for(unsigned int k = 1; k < word_size; k++) {
				if(Grid_str[ind+k] != word[k])
					is_for = false;
				if(Grid_str[ind+k] != word[word_size-k-1])
					is_rev = false;
			}
			if(is_for || is_rev) return true;
		}
		//S
		ind = (loc_i+1-word_size+i)*width+loc_j;
		if(loc_i+1+i >= word_size && loc_i+i < height && ((is_for = (Grid_str[ind] == word[0])) || (is_rev = (Grid_str[ind] == word[word_size-1]))) ) {
			if(is_for) is_rev = Grid_str[ind] == word[word_size-1];
			for(unsigned int k = 1; k < word_size; k++) {
				if(Grid_str[ind+k*width] != word[k])
					is_for = false;
				if(Grid_str[ind+k*width] != word[word_size-k-1])
					is_rev = false;
			}
			if(is_for || is_rev) return true;
		}
		//SE
		ind = (loc_i+1-word_size+i)*width+loc_j+1-word_size+i;
		if(loc_j+1+i >= word_size && loc_j+i < width && loc_i+1+i >= word_size && loc_i+i < height && ((is_for = (Grid_str[ind] == word[0])) || (is_rev = (Grid_str[ind] == word[word_size-1])))) {
			if(is_for) is_rev = Grid_str[ind] == word[word_size-1];
			for(unsigned int k = 1; k < word_size; k++) {
				if(Grid_str[ind+k*(width+1)] != word[k])
					is_for = false;
				if(Grid_str[ind+k*(width+1)] != word[word_size-k-1])
					is_rev = false;
			}
			if(is_for || is_rev) return true;
		}
		//SW
		ind = (loc_i+1-word_size+i)*width+loc_j+word_size-1-i;
		if(loc_j+word_size-1 < width+i && loc_j >= i && loc_i+1+i >= word_size && loc_i+i < height && ((is_for = (Grid_str[ind] == word[0])) || (is_rev = (Grid_str[ind] == word[word_size-1])))) {
			if(is_for) is_rev = Grid_str[ind] == word[word_size-1];
			for(unsigned int k = 1; k < word_size; k++) {
				if(Grid_str[ind+k*(width-1)] != word[k])
					is_for = false;
				if(Grid_str[ind+k*(width-1)] != word[word_size-k-1])
					is_rev = false;
			}
			if(is_for || is_rev) return true;
		}
	}
	return false;
}


// Reverses each line of the grid to flip it horizontally. O(w*h)
void flip_grid_hor(unsigned int width, unsigned int height, std::string& Grid_str) {
	for(unsigned int i = 0; i < height; i++)
		std::reverse(Grid_str.begin()+i*width, Grid_str.begin()+(i+1)*width);
}

// Reverses the entire grid, flipping it horizontally and vertically. O(w*h)
void flip_grid_hor_ver(unsigned int width, unsigned int height, std::string& Grid_str) {
	std::reverse(Grid_str.begin(), Grid_str.end());
}


/** 
 * Swaps letters across the diagonal (think matrix transpose).
 * Is equivalent to rotating grid 90 degrees anticlockwise, then flipping it vertically.
 */
void flip_grid_sqr(unsigned int length, std::string& Grid_str) {
	for(unsigned int i = 0; i < length; i++)
		for(unsigned int j = i+1; j < length; j++)
			std::iter_swap(Grid_str.begin()+i*length+j, Grid_str.begin()+j*length+i);
}

// Returns a boolean indicating whether or not the grid is symmetrical horizontally.
bool is_grid_sym_hor(unsigned int width, unsigned int height, const std::string& Grid_str) {
	for(unsigned int i = 0; i < width; i++)
		for(unsigned int j = 0; j < height/2; j++)
			if(Grid_str[j*width+i] != Grid_str[(height-j-1)*width+i])
				return false;
	return true;
}

// Returns a boolean indicating whether or not the grid is symmetrical vertically.
bool is_grid_sym_ver(unsigned int width, unsigned int height, const std::string& Grid_str) {
	for(unsigned int i = 0; i < height; i++)
		for(unsigned int j = 0; j < width/2; j++)
			if(Grid_str[i*width+j] != Grid_str[i*width+width-j-1])
				return false;
	return true;
}

/**
 * Looks for the first empty ('0') character in the grid and fills it with a valid letter.
 * Then, it checks if the grid has any excluded words in it using is_word_in_check2. If there are
 * no excluded words, it recursively calls itself. Once every letter has been tried, the character
 * is set back to empty and the function is returned, as the recursive calls have already dealt
 * with each of the empty characters following it. If there are no empty characters, then the grid
 * is emplaced into the Grids solution set, as it must be a valid solution. Then, if the grid is
 * not already in the set, its horizontal and/or vertical flips are also emplaced into Grids, and
 * then if the grid is square, a 90 degree rotations and each of its horizontal and/or vertical
 * flips are also emplaced into Grids, as they are also valid solutions.
 */
void check_exclude_fill_emplace(unsigned int width, unsigned int height, std::string& Grid_str,
							    std::unordered_set<std::string>& Grids, Node* exclude_node,
							    bool& find_all, unsigned int i_0) {
	if(!find_all && Grids.size() > 0) return;
	bool is_valid;
	unsigned int grid_size = height*width;
	for(unsigned int i = i_0; i < grid_size; i++) {
		if(Grid_str[i] == '0') {
			for(char c = 'a'; c <= 'z'; c++) {
				Grid_str[i] = c;
				is_valid = true;
				for(Node* exclude_node_it = exclude_node; exclude_node_it; exclude_node_it = exclude_node_it->ptr) {
					if(exclude_node_it->value.find(c) != std::string::npos && is_word_in_check2(width, height, Grid_str, exclude_node_it->value, i)) {
						is_valid = false;
						break;
					}
				}
				if(is_valid) check_exclude_fill_emplace(width, height, Grid_str, Grids, exclude_node, find_all, i);
			}
			Grid_str[i] = '0';
			return;
		}
	}
	unsigned int grids_size = Grids.size();
	Grids.emplace(Grid_str);
	if(!find_all || Grids.size() == grids_size) return;
	flip_grid_hor(width, height, Grid_str);
	Grids.emplace(Grid_str);
	flip_grid_hor_ver(width, height, Grid_str);
	Grids.emplace(Grid_str);
	flip_grid_hor(width, height, Grid_str);
	Grids.emplace(Grid_str);
	flip_grid_hor_ver(width, height, Grid_str);
	if(height != width) return;		
	flip_grid_sqr(width, Grid_str);
	Grids.emplace(Grid_str);
	flip_grid_hor(width, height, Grid_str);
	Grids.emplace(Grid_str);
	flip_grid_hor_ver(width, height, Grid_str);
	Grids.emplace(Grid_str);
	flip_grid_hor(width, height, Grid_str);
	Grids.emplace(Grid_str);
	flip_grid_hor(width, height, Grid_str);
	flip_grid_sqr(width, Grid_str);
}


/**
 * If there are no words left to be included, a check is done for excluded words using
 * is_word_in_check1, then check_exclude_fill_emplace is called.
 * If there are still words to be included, checks are done throughout the grid to determine each
 * possible location that the word could be inserted. If there is symmetry in the grid before a
 * word is inserted, then only one half of the grid is checked, and the solutions that would be
 * missed are accounted for in check_exclude_fill_emplace by putting horizontal and/or vertical
 * flips into the solution set. This leads to many fewer recursive calls, and the cost is minimal.
 * Included words are inserted into the grid horizontally, vertically, and diagonally, and once a
 * word is inserted in a valid location, a recursive call is made using the next Node in the linked
 * list of included words and the updated symmetry of the grid.
 */
void inv_word_search(unsigned int width, unsigned int height,
					 const std::string& Grid_str,
					 std::unordered_set<std::string>& Grids,
					 Node* include_node, Node* exclude_node,
					 bool is_sym_hor, bool is_sym_ver, bool& find_all)
{
	if(!find_all && Grids.size() > 0) return;
	if(!include_node) {
		for(Node* exclude_node_it = exclude_node; exclude_node_it; exclude_node_it = exclude_node_it->ptr)
			if(is_word_in_check1(width, height, Grid_str, exclude_node_it->value)) return;
		std::string Grid_str_copy(Grid_str);
		check_exclude_fill_emplace(width, height, Grid_str_copy, Grids, exclude_node, find_all, 0);
		return;
	}
	bool is_space_for;
	bool is_space_rev;
	bool is_space_for_up;
	bool is_space_rev_up;
	std::string word = include_node->value;
	unsigned int word_size = word.size();
	// Whether or not the word is a palindrome is computed.
	bool is_palindrome = true;
	for(unsigned int i = 0; i < word_size/2; i++)
		if(word[i] != word[word_size-i-1]) 
			is_palindrome = false;
	// Checks that the word is small enough to fit in the grid horizontally.
	if(width >= word_size) {
		unsigned int new_height = height*(2-is_sym_hor)/2+is_sym_hor*height%2;
		unsigned int new_width = (width - word_size+ 1)/(1+(is_sym_ver && is_palindrome))+(is_sym_ver && is_palindrome)*width%2;
		for(unsigned int i = 0; i < new_height; i++) {
			for(unsigned int j = 0; j < new_width; j++) {
				is_space_for = true;
				is_space_rev = (word_size > 1 && !is_sym_hor) || (!is_sym_ver && !is_palindrome);
				for(unsigned int k = 0; k < word_size; k++) {
					if(Grid_str[i*width+j+k] != word[k] && Grid_str[i*width+j+k] != '0')
						is_space_for = false;
					if(is_space_rev && Grid_str[i*width+width-j-k-1] != word[k] && Grid_str[i*width+width-j-k-1] != '0')
						is_space_rev = false;
					if(!(is_space_for || is_space_rev)) break;
				}
				if(is_space_for || is_space_rev) {
					if(is_space_for) {
						std::string Grid_str_copy(Grid_str);
						for(unsigned int k = 0; k < word_size; k++)
							Grid_str_copy[i*width+j+k] = word[k];
						inv_word_search(width, height, Grid_str_copy, Grids, include_node->ptr, exclude_node, (is_sym_hor && i == height/2 && height%2 == 1) || is_grid_sym_hor(width, height, Grid_str_copy), (is_sym_ver && is_palindrome && 2*j == width - word_size) || is_grid_sym_ver(width, height, Grid_str_copy), find_all);
					}
					if(is_space_rev) {
						std::string Grid_str_copy(Grid_str);
						for(unsigned int k = 0; k < word_size; k++)
							Grid_str_copy[i*width+width-j-k-1] = word[k];
						inv_word_search(width, height, Grid_str_copy, Grids, include_node->ptr, exclude_node, is_grid_sym_hor(width, height, Grid_str_copy), is_grid_sym_ver(width, height, Grid_str_copy), find_all);
					}
				}
			}
		}
	}
	/** 
	 * Checks that the word is small enough to fit in the grid vertically, that the word has length
	 * greater than 1 (because inserting a word of length 1 vertically is equivalent to doing so
	 * horizontally), if the grid is square, that it is neither symmetrical horizontally nor
	 * vertically (the 90 degree rotation done in check_exclude_fill_emplace will acocunt for these
	 * solutions), and if only one solution is being found, that no solution has been found yet. 
	 */
	if(height >= word_size && word_size > 1 && (!is_sym_hor || !is_sym_ver || height != width) && (find_all || Grids.size() == 0)) {
		unsigned int new_width = width*(2-is_sym_ver)/2+is_sym_ver*width%2;
		unsigned int new_height = (height - word_size + 1)/(1+(is_sym_hor && is_palindrome))+(is_sym_hor && is_palindrome)*height%2;
		for(unsigned int i = 0; i < new_width; i++) {
			for(unsigned int j = 0; j < new_height; j++) {
				is_space_for = true;
				is_space_rev = !is_sym_ver || (!is_sym_hor && !is_palindrome);
				for(unsigned int k = 0; k < word_size; k++) {
					if(Grid_str[(j+k)*width+i] != word[k] && Grid_str[(j+k)*width+i] != '0')
						is_space_for = false;
					if(is_space_rev && Grid_str[(height-j-k-1)*width+i] != word[k] && Grid_str[(height-j-k-1)*width+i] != '0')
						is_space_rev = false;
					if(!(is_space_for || is_space_rev)) break;
				}
				if(is_space_for || is_space_rev) {
					std::string Grid_str_copy(Grid_str);
					if(is_space_for) {
						for(unsigned int k = 0; k < word_size; k++)
							Grid_str_copy[(j+k)*width+i] = word[k];
						inv_word_search(width, height, Grid_str_copy, Grids, include_node->ptr, exclude_node, (is_sym_hor && is_palindrome && 2*j == height - word_size) || is_grid_sym_hor(width, height, Grid_str_copy), (is_sym_ver && i == width/2 && width%2 == 1) || is_grid_sym_ver(width, height, Grid_str_copy), find_all);
					}
					if(is_space_rev) {
						for(unsigned int k = 0; k < word_size; k++)
							Grid_str_copy[(height-j-k-1)*width+i] = word[k];
						inv_word_search(width, height, Grid_str_copy, Grids, include_node->ptr, exclude_node, is_grid_sym_hor(width, height, Grid_str_copy), is_grid_sym_ver(width, height, Grid_str_copy), find_all);
					}
				}
			}
		}
	}
	// Checks that the word is small enough to be put into the grid diagonally, that the word has
	// length greater than 1, and if only one solution, that none have been found yet.
	if(height >= word_size && width >= word_size && word_size > 1 && (find_all || Grids.size() == 0)) {
		for(unsigned int i = 0; i < height - word_size + 1; i++) {
			for(unsigned int j = 0; j < width - word_size + 1; j++) {
				is_space_for = true;
				is_space_rev = !is_sym_ver || (!is_sym_hor && !is_palindrome);
				is_space_for_up = !is_sym_hor || (!is_sym_ver && !is_palindrome);
				is_space_rev_up = !is_palindrome || is_space_rev || is_space_for_up;
				for(unsigned int k = 0; k < word_size; k++) {
					if(Grid_str[(i+k)*width+j+k] != word[k] && Grid_str[(i+k)*width+j+k] != '0')
						is_space_for = false;
					if(is_space_rev && Grid_str[(i+k)*width+width-j-k-1] != word[k] && Grid_str[(i+k)*width+width-j-k-1] != '0')
						is_space_rev = false;
					if(is_space_for_up && Grid_str[(height-i-k-1)*width+j+k] != word[k] && Grid_str[(height-i-k-1)*width+j+k] != '0')
						is_space_for_up = false;
					if(is_space_rev_up && Grid_str[(height-i-k-1)*width+width-j-k-1] != word[k] && Grid_str[(height-i-k-1)*width+width-j-k-1] != '0')
						is_space_rev_up = false;
					if(!(is_space_for || is_space_rev || is_space_for_up || is_space_rev_up)) break;
				}
				if(is_space_for || is_space_rev || is_space_for_up || is_space_rev_up) {
					std::string Grid_str_copy(Grid_str);
					if(is_space_for) {
						for(unsigned int k = 0; k < word_size; k++)
							Grid_str_copy[(i+k)*width+j+k] = word[k];
						inv_word_search(width, height, Grid_str_copy, Grids, include_node->ptr, exclude_node, is_grid_sym_hor(width, height, Grid_str_copy), is_grid_sym_ver(width, height, Grid_str_copy), find_all);
					}
					if(is_space_rev) {
						for(unsigned int k = 0; k < word_size; k++)
							Grid_str_copy[(i+k)*width+width-j-k-1] = word[k];
						inv_word_search(width, height, Grid_str_copy, Grids, include_node->ptr, exclude_node, is_grid_sym_hor(width, height, Grid_str_copy), is_grid_sym_ver(width, height, Grid_str_copy), find_all);
					}
					if(is_space_for_up) {
						for(unsigned int k = 0; k < word_size; k++)
							Grid_str_copy[(height-i-k-1)*width+j+k] = word[k];
						inv_word_search(width, height, Grid_str_copy, Grids, include_node->ptr, exclude_node, is_grid_sym_hor(width, height, Grid_str_copy), is_grid_sym_ver(width, height, Grid_str_copy), find_all);
					}
					if(is_space_rev_up) {
						for(unsigned int k = 0; k < word_size; k++)
							Grid_str_copy[(height-i-k-1)*width+width-j-k-1] = word[k];
						inv_word_search(width, height, Grid_str_copy, Grids, include_node->ptr, exclude_node, is_grid_sym_hor(width, height, Grid_str_copy), is_grid_sym_ver(width, height, Grid_str_copy), find_all);
					}
				}
			}
		}
	}
	return;
}


int main(int argc, char* argv[]) {
	/**
	 * Declares variables while it checks for the correct number of command line arguments,
	 * verifies that the given input file exists, and confirms the command line argument
	 * formatting. 
	 */
	if(argc < 4) {
		std::cerr << "ERROR 1: NOT ENOUGH CMD LINE ARGS." << std::endl << "FORMAT AS:   input_file.txt output_file.txt solution_type" << std::endl;
		return 1;
	}
	std::ifstream file_in(argv[1]);
	if(!file_in) {
		std::cerr << "ERROR 2: INPUT FILE DOES NOT EXIST." << std::endl;
		return 2;
	}
	std::ofstream file_out(argv[2]);
	std::string s_arg(argv[3]);
	bool find_all;
	if(s_arg == "find_all_solutions")
		find_all = true;
	else if(s_arg == "one_solution")
		find_all = false;
	else {
		std::cerr << "ERROR 3: INCORRECT solution_type GIVEN." << std::endl << "MUST BE EITHER \"find_all_solutions\" OR \"one_solution\"" << std::endl;
		return 3;
	}
	std::string cmd = "";
	unsigned int width = 0;
	unsigned int height = 0;
	Node* include_node = NULL;
	Node* include_node_head = NULL;
	Node* exclude_node = NULL;
	Node* exclude_node_head = NULL;
	Node* temp_node = NULL;
	// unordered_set will hold the solutions: no duplicate checking necessary.
	std::unordered_set<std::string> Grids;
	std::string included_letters = "";


	// Grid width and height input and their validity in confirmed.
	file_in >> width;
	file_in >> height;
	if(width == 0 || height == 0) {
		std::cerr << "ERROR 4: INVALID GRID DIMENSIONS IN INPUT FILE." << std::endl << "FIRST LINE MUST BE" << std::endl << "width height" << std::endl << "WHERE width AND height ARE INTS GREATER THAN 0." << std::endl;
		return 4;
	}
	
	/**
	 * while loop reads in words from inout file. Included words are put into a linked list sorted
	 * from longest to shortest, as it is more efficient to attempt placing the larger words first
	 * because this leaves less open spaces and thus less attempted starting points for words that
	 * follow it, leading to fewer recursive calls of inv_word_search. Excluded words are put into
	 * a linked list sorted from shortest to longest, as the smaller words are not only more likely
	 * to appear, but also take less time to check, leading to less time spent checking excluded
	 * words. Also, format checking for included/excluded words is done.
	 */
	while(file_in >> cmd) {
		if(cmd == "+") {
			file_in >> cmd;
			// Checks that the input string is the proper format and adds any new letters in the
			// word to the included_letters string.
			for(unsigned int i = 0; i < cmd.size(); i++) {
				if(included_letters.find(cmd[i]) == std::string::npos) {
					included_letters += cmd[i];
					if(cmd[i] < 'a' || cmd[i] > 'z') {
						std::cerr << "ERROR 6: ONLY LOWERCASE LETTERS CAN BE IN INCLUDED WORDS." << std::endl << cmd << " CONTAINS " << cmd[i] << ", AN INAVLID CHARACTER." << std::endl;
						del_linked_list(include_node_head);
						del_linked_list(exclude_node_head);
						return 6;
					}
				}
			}
			if(!include_node_head || cmd.size() > include_node_head->value.size()) {
				include_node = new Node;
				include_node->value = cmd;
				include_node->ptr = include_node_head;
				include_node_head = include_node;
			}
			else {
				for(include_node = include_node_head; include_node; include_node = include_node->ptr) {
					if(!include_node->ptr || cmd.size() >= include_node->ptr->value.size()) {
						temp_node = include_node->ptr;
						include_node->ptr = new Node;
						include_node = include_node->ptr;
						include_node->value = cmd;
						include_node->ptr = temp_node;
						break;
					}
				}
			}
		}
		else if(cmd == "-") {
			file_in >> cmd;
			for(unsigned int i = 0; i < cmd.size(); i++) {
				if(cmd[i] < 'a' || cmd[i] > 'z') {
					std::cerr << "ERROR 7: ONLY LOWERCASE LETTERS CAN BE IN EXCLUDED WORDS." << std::endl << cmd << " CONTAINS " << cmd[i] << ", AN INAVLID CHARACTER." << std::endl;
					del_linked_list(include_node_head);
					del_linked_list(exclude_node_head);
					return 7;
				}
			}
			if(cmd.size() <= width || cmd.size() <= height) {
				if(!exclude_node_head || cmd.size() <= exclude_node_head->value.size()) {
					exclude_node = new Node;
					exclude_node->value = cmd;
					exclude_node->ptr = exclude_node_head;
					exclude_node_head = exclude_node;
				}
				else {
					for(exclude_node = exclude_node_head; exclude_node; exclude_node = exclude_node->ptr) {
						if(!exclude_node->ptr || cmd.size() < exclude_node->ptr->value.size()) {
							temp_node = exclude_node->ptr;
							exclude_node->ptr = new Node;
							exclude_node = exclude_node->ptr;
							exclude_node->value = cmd;
							exclude_node->ptr = temp_node;
							break;
						}
					}
				}
			}
		}
		else {
			std::cerr << "ERROR 5: MISSING +/- IN INPUT FILE." << std::endl << "MUST HAVE + BEFORE EACH WORD TO BE INCLUDED, - BEFORE THOSE TO BE EXCLUDED." << std::endl;
			del_linked_list(include_node_head);
			del_linked_list(exclude_node_head);
			return 5;
		}
	}

 	
 	// Creates an "empty" (filled with 0s) grid to pass to the recursive function, but
 	// first checks that there is enough area in the grid for every included letter.
	std::string Grid_str(width*height, '0');
	if(width*height >= included_letters.size())
		inv_word_search(width, height, Grid_str, Grids, include_node_head, exclude_node_head, true, true, find_all);

	// Prints out the results from the Grids unordered_set if there are any.
	if(Grids.size() > 0) {
		if(find_all)
			file_out << Grids.size() << " Solution(s)" << std::endl;
		for(typename std::unordered_set<std::string>::iterator it = Grids.begin(); it != Grids.end(); ++it) {
			file_out << "Board: " << std::endl;
			for(unsigned int i = 0; i < height; ++i) {
				file_out << "  ";
				for(unsigned int j = 0; j < width; ++j)
					file_out << (*it)[i*width+j];
				file_out << std::endl;
			}
		}
	}
	else file_out << "No solutions found" << std::endl;

	del_linked_list(include_node_head);
	del_linked_list(exclude_node_head);
	
	return 0;
}