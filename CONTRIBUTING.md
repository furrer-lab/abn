# Contributing to abn

This outlines how to propose a change to abn.

## Fixing typos

If you spot typos or other language related errors in the documentation you may directly us the GitHub web interface to edit the source file.

## Beyond typos

If you would like to contribute bigger changes please always refer to the
[issue board](https://github.com/furrer-lab/abn/issues).

In case you simply consider contributing (that's awesome!) head over to the
[issue board](https://github.com/furrer-lab/abn/issues) where you will find
open issues  with the label
https://github.com/furrer-lab/abn/labels/help%20wanted.
Feel free to comment on them if it is not clear to you what exactly the issue
is about.
If you have found an open issue that you would like to fix simply follow the
[Pull request process](#pull-request-proces) that is described below.

In case you already have some specific contribution you would like to add,
please first check (also in the already closed issues) if you can find existing
issues on [issue board]([issue board](https://github.com/furrer-lab/abn/issues)
that are related to the contribution you plan to make.
If yes, please leave a comment on the related issue.
If you cannot find an issue related to the contribution you plan to make, please
create a [new issue](https://github.com/furrer-lab/abn/issues/new).
For some inspiration on how to write a great issue refer to the
[tidyteam code review prinicples](https://code-review.tidyverse.org/issues/).

Once your foreseen contribution is documented in an issue, feel free to start
the [Pull request process](#pull-request-proces).


### Pull request process

*   Fork the package and clone onto your computer.
    You might use `usethis::create_from_github("furrer-lab/abn", fork = TRUE)`
    to do so.
*   `abn` ships with a development environment provided by `renv`...
