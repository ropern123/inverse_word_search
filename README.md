HOMEWORK 6 CONTEST: INVERSE WORD SEARCH


NAME:  < Nicholas Roper >


OPTIMIZATIONS:
I sort the list of words that need to be included by largest to smallest, as the
larger words take up more space and leave less places where the smaller words might
go.

I sort the list of words that need to be excluded by smallest to largest, as the
smaller words are faster to check and more likely to be in the grid.

My grid is a single string of size width\*height that the program treats as a
2d grid, i.e. I iterate through it using Grid_str[i\*width+j], with i being the
row and j being the column.

I first try to add words horizontally, then vertically, then diagonally. In all of these,
I also check to see if the words can be added in reverse at the same time.

I check for symmetry to reduce the number of places words are put to reduce the number
of recursive calls, as many of the possible solutions are usually the same as some of
the others, but they have been flipped horizontally or vertically or both.

As an addition to symmetry, I check if words being put into the grid are palidromes so that
they do not need to be added both forwards and backwards.

I check to see if the grid is a square grid, as there are more symmetries in this case,
as each solution can be the same as many of the others, except tilted 90 degrees
(i.e. ab  -->  90 deg   -->  bd )
(     cd  -->  becomes  -->  ac ).

I check to see of the words that are being added are 1 letter long, as if they are, they
do not need to be added horizontally, vertically, and diagonally, as they will be the same
for each, so they are instead added only horizontally.

My first function that checks to see if a word is in the grid (is_word_in_check1) only checks
to the East, Southeast, South, and Southwest, as it also checks in reverse at the same time,
and so doesn't have to check through as much of the grid.

My second function that checks to see if a word is in the grid (is_word_in_check2) is only
called after a single character has been added, so it only checks the positions of the right
length that might have become an excluded word with the addition of that letter, greatly
reducing the amount of the board that is checked. It checks in the same directions as
is_word_in_check1 and also both forwards and reverse to get every possible word location. 

I used unordered sets to hold my grids, as they cannot hold duplicates (so no checking needs
to be done) and they have average O(1) emplace.





TEST CASE SUMMARY:
It runs all of the 8 given test cases in under 0.06s.

Test case 1: ropern_1.txt (find_all_solutions) Creates a 7x7 grid that is symmetric horizontally and vertically with only two solutions. It has 8 included words of length 7 that fill up each quadrant of the grid except the lower right 3x3, 60 excluded
words of length 2 and 10 excluded words of length 1, all of which are necessary to bring the number of solutions down to only
two. My program runs it in (on cygwin, optimized with -O3) .723s.

Test case 2: ropern_2.txt (one_solution) Creates a large grid with every letter excluded except for n (as it is in the middle and will take a long time if the excluded words are sorted forwards or backwards). If the program does not check for excluded
words at each step in the process of filling or only use n for its alphabet, it will not finish as there are too many possible
combinations that will be tried. My program run it in .06s.