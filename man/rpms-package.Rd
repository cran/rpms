% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rpms.R
\docType{package}
\name{rpms-package}
\alias{rpms-package}
\title{Recursive Partitioning for Modeling Survey Data (rpms)}
\description{
This package provides a function \code{rpms} to produce an \code{rpms} object 
and method functions that operate on them. 
The \code{rpms} object is a representation of a regression tree achieved
by recursively partitioning the dataset, fitting the specified linear model
on each node separately.
The recursive partitioning algorithm has an unbiased variable selection
and accounts for the sample design.
The algorithm accounts for one-stage of stratification and clustering as
well as unequal probability of selection.
There are also functions for producing random forest estimator 
(a list of \code{rpms} objects), a boosted regression tree and tree 
based zero-inflated model.
}
