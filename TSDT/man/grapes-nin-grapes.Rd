\name{\%nin\%}
\alias{\%nin\%}
\title{\%nin\%}

\usage{a \%nin\% b}

\arguments{
  \item{a}{Any R object for which the binary operator \%in\% is
    defined. This would include many built-in R primitives.}
  \item{b}{Any R object for which the binary operator \%in\% is
    defined. This would include many built-in R primitives.}
}

\description{Negation of the built-in \%in\% operator. \%nin\% is
  a short-hand for !( a \%in\% b ).
}

\seealso{\link{\%in\%}}
  
\examples{
# 4 is not an element in {5,6,7}.
4 \%nin\% 5:7  # Evaluates to TRUE

# 4 is an element in {4,5,6,7}.
4 \%nin\% 4:7  # Evaluates to FALSE
}
