# Contributing to `abn`

This outlines how to propose a change to `abn`.

## The "Pull request" process

In general all edits to the `abn` package go through the process of a
[Pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request).

### Fixing typos

If you spot typos or other language related errors in the documentation you may
directly use the GitHub web interface to edit the source file.
Simply click on the edit :pen: button and Github will guide you through the
process of forking the repository and creating a pull request.

### Beyond typos

If you would like to contribute bigger changes please always refer to the
[issue board](https://github.com/furrer-lab/abn/issues) first!

> [!NOTE]
> With the exception of simple typos, never open a pull request without an
> existing issue that it relates to!
>
> If you want to link a pull request to an issue simply copy the link to the
> issue into the description of the pull request.
> If you are unsure how to to that have a look at the
> official [Github documentation for pull requests](https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue) for further details.

#### If you want to help,

but have nothing specific in mind, simply head over to the [issue board](https://github.com/furrer-lab/abn/issues) where
you will find open issues labelled with
https://github.com/furrer-lab/abn/labels/help%20wanted.
Feel free to comment on them if it is not clear to you what exactly the issue is about.
Issues with the label
https://github.com/furrer-lab/abn/labels/good%20first%20issue might be good
starting points if you are new to `abn`.

#### If you have a specific contribution
that you would like to add first check the
[issue board](https://github.com/furrer-lab/abn/issues)
(also see the
[closed issues](https://github.com/furrer-lab/abn/issues?q=is%3Aissue+is%3Aclosed))
if you can find existing issues related to the contribution you plan to make.
If yes, leave a comment on the related issue.
If no, create a [new issue](https://github.com/furrer-lab/abn/issues/new) describing
the feature you want to implement.
For some inspiration on how to write a great issue refer to the
[tidyteam code review prinicples](https://code-review.tidyverse.org/issues/).

Once the foreseen edits are documented in an issue, create a pull request
(to see how simply follow
[Github's instructions on creating a Pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request))
and add your changes to the related branch.

## Development

We use [renv](https://rstudio.github.io/renv/index.html) for our development
environment and would encourage you to do the same.
Therefore, `abn` directly provides a development environment and we suggest the
following workflow for development:

*   Fork the package and clone it onto your computer.
    
    _Note:_ You might do this with `usethis::create_from_github("furrer-lab/abn", fork = TRUE)`.

*   Install `renv` with `install.packages("renv")`
*   Restore the development environment with `renv::restore()`
*   Create a new branch for your pull request.
    
    _Note:_ You might use `usethis::pr_init("explicative-branch-name")`.
    
    _See [Testing](#testing) for flags you can set in the branch name_

*   Edit the code and commit your changes.
  
    _Note:_ With `usethis` you would run: `usethis::pr_push()` and the follow 
    the prompts.
    The title of your pull request should briefly describe the change.
    The body of your pull request should contain `Closes #issue-number`, where
    `#issue-number` is the identifier of the related issue.

## Testing

The `abn` package uses several testing pipelines.
When developing (i.e. not on the `main` branch) we recommend running the
[quick-testthat.yml](https://github.com/furrer-lab/abn/blob/24-documentation-of-the-testing-procedure-noT/.github/workflows/quick-testthat.yml) action.
This pipeline runs by default on push event on all branches other than `main`
and executes all the test defined withing the `abn` package by calling
[`devtools::test`](https://devtools.r-lib.org/reference/test.html).
It does so in a docker container from the
[r-containers](https://github.com/furrer-lab/r-containers) package which builds
container images that are setup with all the required dependencies.

By default, a container based on `Debian` using `clang` as compiler front-end
and the current R development version version (`debian/clang/devel`) is used.
However, [other configurations are available](https://github.com/orgs/furrer-lab/packages?repo_name=r-containers).

If you want to change the container image to use for testing your branch against,
you can append the desired configuration to the name of your branch with the
prefix `'__'`.
Hence, a branch named `my-awesome-branch__debian/gcc/release` would use the
container image based on Debian, using `gcc` and the current release version of
R to run the tests.

You can also completely opt-out of testing by adding the flag `noT` to your
branch name.
It does not matter where you put the string `noT` in the branch name, so a
branch named `"monoTestingBranch"` will also make you opt-out of testing.

_Note:_ If you start a commit message with `"noT"` then the testing pipeline
will also not run for this commit. In a commit message, however, the message
has to start with `noT`, it will be ignored otherwise.

### Memory Checks

By default our testing pipelines do not perform checks for proper memory usage.
However, if you are working on the C code, you have to option to run computationally
more expensive (i.e. slower) memory checks.

To run such checks, your development branch needs to have an open Pull request and you
simply have to label the request with https://github.com/furrer-lab/abn/labels/memory%3A%3Acheck
