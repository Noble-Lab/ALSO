
#  Do this: https://stackoverflow.com/questions/12226846/count-letter-differences-of-two-strings
#+ - https://www.infoworld.com/article/3340120/how-to-run-python-in-r.html
#+ - https://www.w3schools.com/python/ref_func_zip.asp
#+ - https://www.w3schools.com/python/ref_func_dict.asp
#+ - https://stackoverflow.com/questions/7818970/is-there-a-dictionary-functionality-in-r
#+ - https://stackoverflow.com/questions/9281323/zip-or-enumerate-in-r
#+ - https://stackoverflow.com/questions/159720/what-is-the-naming-convention-in-python-for-variable-and-function-names
#+ - https://www.infoworld.com/article/3340120/how-to-run-python-in-r.html


library(magrittr)
library(reticulate)

use_condaenv(
    "/Users/kalavattam/miniconda3/envs/r41_env",
    required = TRUE
)
system("which python")
system("python --version")

string_1 <- "IGBDKYFHARGNYDAA"
string_2 <- "KGADKYFHARGNYEAA"
# string_1 <- "Y2xpZW50X2lkOjRjNDY1OGM1NzM"
# string_2 <- "Y2xpZW50X2lkOjRjNOY1OGM1NzM"
# string_1 <- 'bbc'
# string_2 <- 'abd'
# string_1 <- 'abcji'
# string_2 <- 'abdjk'

python_code <- 
"
def compare(string_1, string_2, no_match_c=' ', match_c='|'):
    if len(string_2) < len(string_1):
        string_1, string_2 = string_2, string_1
    result = ''
    n_diff = 0
    for c1, c2 in zip(string_1, string_2):
        if c1 == c2:
            result += match_c
        else:
            result += no_match_c
            n_diff += 1
    delta = len(string_2) - len(string_1)
    result += delta * no_match_c
    n_diff += delta
    return (result, n_diff)


def main():
    string_1 = r.string_1
    string_2 = r.string_2
    result, n_diff = compare(string_1, string_2, no_match_c='*')
    print(string_1)
    print(result)
    print(string_2)
    print('''%d difference(s).\n''' % n_diff)
    

main()
"

python_code %>% reticulate::py_run_string()
