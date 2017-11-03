## Test environments
* local OS X install, R 3.4.2
* ubuntu 12.04 (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---


* In reply to the first round comment on the first submission  

> Please only use + file LICENSE for restrictions to the GPL-3 
and only include the restrictions.

We need add no restrictions to GPL-3 so The **License** field in `DESCRIPTION` 
was modified to `GPL-3` and the `LICENSE` file was removed. 
In addition, a file `cran-comments.md` was added to track comments.  

* Second round

> Thanks, please elaborate in your description the used methods (e.g. which tests).   
  Please add a reference in the 'Description' field of your DESCRIPTION file reference in the form  
  authors (year) <doi:...>  
  authors (year) <arXiv:...>  
  authors (year, ISBN:...)  
  with no space after 'doi:', 'arXiv:' and angle brackets for auto-linking.  
  Please explain the acronym PCR in your description.   
  Please fix and resubmit.  

We rephrased the **Description** field of the `DESCRIPTION` file to elaborate 
on the methods, explain the acronym PCR and add two references  

* v1.0.1 submission

Correcting a typo in the DESCRIPTION file
