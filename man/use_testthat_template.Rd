% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_template.R
\name{use_testthat_template}
\alias{use_testthat_template}
\title{Create tests/testthat/test-*.R test file}
\usage{
use_testthat_template(name)
}
\arguments{
\item{name}{A character string providing the name of the R function that you
want to test. The name will be used to create the file name, i.e.,
\verb{tests/testthat/test-\{name\}.R}. If \code{name} is not specified, the function
will not be able to create a file name and an error will be returned.}
}
\value{
If successful, this function invisibly returns \code{TRUE} to allow for the
chaining of commands. If the function is unsuccessful, an error message is
returned.

Two messages are also returned from the usethis package, which is used by
this function. The first states where the FIMS project is on your computer.
The second states the file path of the newly created file. The file will not
be automatically opened.
}
\description{
This helper function creates a new test file for testthat using a template
available within the templates folder of this package. The test file is
created under the "tests/testthat" directory with a name based on the input
argument \code{name}. If the test file already exists, an error message is
returned and no changes are made.
}
\details{
There are three minimum testing criteria for FIMS, which should be validated
for every R function within the package. The template file sets up a
section for each of the three following test criteria:
\enumerate{
\item Input and output correctness (IO correctness): ensure that the function
behaves as expected with correct inputs and returns the expected outputs.
\item Edge-case handling (Edge handling): validate the function's performance
with invalid inputs and unusual scenarios.
\item Built-in errors and warnings (Error handling): confirm that appropriate
error and warning messages are triggered under exceptional conditions.
}

Above every expectation within the test file there should be a call to
\verb{#' @description} that fits on one line describing the test. The information
will be used in the bookdown report of the testing results.
}
\examples{
# Create a new test file named "test-new_function" for `new_function()`
\dontrun{
FIMS:::use_testthat_template("new_function")
}
}
