# Coding Recommendations

There are a few things which, although far from guaranteeing a swift merge of your pull request, _might_ well smooth its passage.


## Coding Practice
 - check system calls for errors
 - use ```kstring``` utilities for string manipulations rather than homebrewing or adding another dependancy
 - avoid using yoda expressions to avoid accidental assignment (e.g. avoid ```if ( some_func() == some_var)``` and ```if ( 0 == some_var)``` rather than  ```if (some_var == some_func())``` and ```if (some_var == 0)```) as we have warnings turned on to help prevent that 
 - avoid introducing regular expressions (even though they get functionality quickly they often cause problems later)
 - use ```const``` appropriately for read-only data, especially in an API

##Style
Fitting in with the local code is no bad thing.

If in doubt or if the local style is inconsistent:
 - no tabs
 - four space indentation
 - no trailing space
 - braces - in function definitions branches should be matching in the same column i.e. ```{``` on line after its name
 - try not to go over line lengths of 132

You'll find plenty of counter examples to the above suggestions (and there's probably not currently much appetite for a lot of whitespace churn).

##Git & GitHub
 - rebasing pull requests to the latest state of develop to make merging easier is appreciated _unless_ the PR is being actively considered and worked on.
