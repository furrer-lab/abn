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

The `abn` package uses a comprehensive multi-layered testing infrastructure designed to ensure code quality while being efficient during development. Understanding when and how different tests run is essential for effective contribution.

### Automatic Testing During Development

#### Development Tests (`development_run.yml`)
When you create a pull request to any branch (except `main` and release branches), the package automatically runs **intelligent testing**:

- **Smart test selection**: Only runs tests related to the files you modified
  - If you modify files in `src/`: runs all tests (C code changes can affect anything)
  - If you modify files in `R/`: runs only tests for the specific R files changed
  - If you modify other files: may skip testing entirely

- **Default environment**: Tests run in a Docker container with Debian, `clang` compiler, and R development version (`debian/clang/devel`)

- **Custom environments**: You can test against different configurations by appending `__<config>` to your branch name:
  ```
  my-feature-branch__debian/gcc/release    # Use Debian + GCC + R release
  my-feature-branch__fedora/clang/devel    # Use Fedora + clang + R devel
  ```
  [See available configurations](https://github.com/orgs/furrer-lab/packages?repo_name=r-containers).

- **Opt-out option**: Skip testing by including `noT` anywhere in your branch name or starting a commit message with `noT`.

### On-Demand Testing via Labels

For thorough testing beyond the basic development tests, apply labels to your pull request to trigger specialized checks:

| Label | Purpose | When to Use |
|-------|---------|-------------|
| `CRAN::check` | Full R CMD check across multiple OS/compiler combinations | Before merging significant changes, especially for CRAN submission |
| `memory::check` | Memory leak detection with Valgrind | When modifying C/C++ code in `src/` |
| `linting::check` | Code style and linting validation | For code cleanup or style improvements |
| `URL::check` | Validates all URLs in documentation | When updating documentation or links |
| `pkgdown::check` | Builds package website | When modifying documentation or vignettes |
| `CompVignettes::build` | Rebuilds computational vignettes | When vignette code or examples change |

**How to trigger**: Simply add the appropriate label to your pull request. The workflow will automatically run and update the label to indicate success (`::passed`) or failure (`::failed`).

### Release Testing

#### Release Branch Testing (`release_run.yml`)
When working on release branches (named `release-X.Y.Z`):
- Automatic linting checks
- Preparation for CRAN submission

#### Installation Testing
On pushes to `main`, the package is automatically tested for installation across:
- Ubuntu (latest)
- Fedora (latest) 
- macOS (setup workflows)
- Windows (setup workflows)

### Testing Strategy by Contribution Type

#### Small Bug Fixes or Documentation Changes
- Development tests are usually sufficient
- Consider `linting::check` for style consistency

#### New Features or Significant Changes
1. Start with development tests during development
2. Before requesting review, add `CRAN::check` label
3. If modifying C code, add `memory::check` label
4. If changing vignettes, add `CompVignettes::build` label

#### Pre-Release Preparation
- All on-demand tests should pass
- Full CRAN check across all platforms
- Memory checks for any C code changes
- Documentation and vignette rebuilds

## Version Bump and Release

The `abn` package follows a structured release process using **Semantic Versioning (SemVer)** and automated workflows. Understanding this process is essential for maintainers and contributors involved in preparing releases.

### Semantic Versioning

The `abn` package uses semantic versioning with the format `MAJOR.MINOR.PATCH`:

- **MAJOR** (`X.0.0`): Breaking changes that are not backward compatible
- **MINOR** (`X.Y.0`): New functionality added in a backward compatible manner  
- **PATCH** (`X.Y.Z`): Backward compatible bug fixes and minor improvements

For example, the package is currently at version `3.1.10`, indicating it's in the third major version with ongoing minor updates and patches.

### Release Process Overview

The release process is fully automated through GitHub Actions and consists of three main phases:

1. **Initiation**: Triggered by creating a release candidate tag
2. **Testing & Validation**: Comprehensive testing of the release candidate
3. **Publishing**: Creating the final release and distributing artifacts

### Step-by-Step Release Process

#### 1. Prepare for Release

Before initiating a release:

- Ensure all planned features/fixes are merged into `main`
- Verify that all tests pass on `main` branch
- Update documentation and vignettes if needed
- Review and update `NEWS.md` manually if necessary (though it's auto-generated)

#### 2. Initiate Release Process

**Trigger**: Create a release candidate tag with the format `X.Y.Z-rc`

```bash
# Example for version 3.1.11
git tag 3.1.11-rc
git push origin 3.1.11-rc
```

**What happens automatically** (`initiate_version_release.yml`):

1. **Version Detection**: Extracts the target version (e.g., `3.1.11` from `3.1.11-rc`)
2. **Changelog Generation**: Updates `NEWS.md` using git-chglog with changes since the last release
3. **File Updates**:
   - `DESCRIPTION`: Updates version and date
   - `configure.ac`: Updates version number
   - `README.md`: Updates status badge branch references
4. **Release Branch Creation**: Creates a `release-3.1.11` branch
5. **Pull Request Creation**: Opens a draft PR from the release branch to `main`
6. **Automatic Labeling**: Adds testing labels (`CRAN::check`, `memory::check`, `linting::check`, `URL::check`)

#### 3. Release Validation Phase

Once the release PR is created, **all comprehensive tests run automatically**:

- **CRAN Checks**: Full `R CMD check` across multiple OS/compiler combinations
- **Memory Checks**: Valgrind testing for memory leaks (if C code changed)
- **Linting**: Code style validation
- **URL Validation**: Checks all documentation links
- **Release Branch Testing**: Basic linting on the release branch

**Manual Steps Required**:
1. Review the automatically generated changelog in `NEWS.md`
2. Make any necessary manual edits to the release PR
3. Ensure all automated tests pass (labels show `::passed`)
4. Convert the draft PR to ready for review
5. Get the release approved by maintainers

#### 4. Publishing the Release

**Option A: Merge Release PR** (Recommended)

When the release PR is merged, `publish_version_merging.yml` automatically:
1. Builds the R package tarball
2. Generates pkgdown documentation
3. Creates a GitHub release with artifacts
4. Deploys documentation to GitHub Pages
5. Labels the PR as `release::published`

**Option B: Push Final Tag** 

Alternatively, push the final version tag directly:

```bash
git tag 3.1.11
git push origin 3.1.11
```

This triggers `publish_version_enforce.yml` with the same publishing actions.

### Making Minor Changes During Release

#### Before the Release PR is Merged

**Small fixes** (typos, documentation):
1. Make changes directly to the release branch (`release-X.Y.Z`)
2. Push to the release branch
3. Changes will be included when the PR is merged

**Significant changes**:
1. Consider if the version number needs adjustment
2. Update the release branch accordingly
3. Rerun testing by adding/removing labels to the PR

#### After Release is Published

**Patch releases** (urgent bug fixes):
1. Create a new patch version (e.g., `3.1.11` → `3.1.12`)
2. Follow the same release process

**Post-release minor changes**:
1. Make changes on `main` branch
2. Include in next scheduled release

### CRAN Submission Process

After a successful GitHub release:

1. **Download the Package**: Get the built `.tar.gz` from the GitHub release
2. **Final CRAN Check**: Run additional CRAN checks locally:
   ```bash
   R CMD check --as-cran abn_X.Y.Z.tar.gz
   ```
3. **Submit to CRAN**: Upload through the [CRAN submission portal](https://cran.r-project.org/submit.html)
4. **Monitor Submission**: Check for automated CRAN feedback and address any issues

### Best Practices

- **Release Frequency**: Plan releases according to accumulated changes and user needs
- **Version Planning**: Use patch versions for bug fixes, minor versions for new features
- **Testing**: Always let all automated tests complete before finalizing a release  
- **Documentation**: Ensure pkgdown site builds successfully before release
- **Communication**: Announce significant releases to users through appropriate channels

### Troubleshooting Releases

- **Failed Tests**: Address test failures before proceeding with release
- **Build Issues**: Check container configurations and dependencies
- **Version Conflicts**: Ensure version numbers are consistent across all files
- **CRAN Rejections**: Address CRAN feedback and create a patch release if needed
