\name{lna_dataset_declaration}
\alias{lna_dataset_declaration}
\alias{validate_lna_dataset_declaration}
\alias{write_lna_dataset_declaration}
\alias{read_lna_dataset_declaration}
\title{Create and read the Neuroarchive dataset declaration}
\usage{
lna_dataset_declaration(plan, participants_receipt = NULL)

validate_lna_dataset_declaration(declaration, plan_hash = NULL)

write_lna_dataset_declaration(declaration, path, overwrite = FALSE)

read_lna_dataset_declaration(path)
}
\arguments{
\item{plan}{A validated LNA dataset plan.}

\item{participants_receipt}{Optional receipt for an explicitly selected
participant row and column subset.}

\item{declaration}{A declaration object or decoded JSON object.}

\item{plan_hash}{Optional expected plan hash.}

\item{path}{Archive root or a path named \code{lna_dataset.json}.}

\item{overwrite}{Whether an existing declaration may be replaced.}
}
\value{
The constructor, validator, and reader return a validated JSON-safe declaration.
The writer returns its path invisibly.
}
\description{
\code{lna_dataset.json} declares the Neuroarchive exchange profile, manifest
authority, metadata and privacy policy, and the \code{.lna.h5} extension and
media type.
}
\details{
The declaration states that the archive uses a BIDS-shaped envelope while
\code{.lna.h5} is a Neuroarchive extension not recognized by the BIDS standard.
This makes the divergence discoverable without filename sniffing.
}
