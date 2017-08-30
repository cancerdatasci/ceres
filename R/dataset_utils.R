#' remove rows from a matrix that include any NAs
#'
#' @param x a matrix
#' @param indices logical, return indices to remove instead
#' @return the matrix with rows removed, unless indices is TRUE
#' @export
remove.rows.all.nas <- function(x, indices=FALSE) {
    i <- plyr::aaply(x, 1, function(r) {all(is.na(r))})
    if (indices) {
        return(which(i))
    } else {
        return(x[ ! i, ])
    }
}

#' remove rows from a matrix that are entirely NAs
#'
#' @inheritParams remove.rows.all.nas
#' @return the matrix with rows removed, unless indices is TRUE
#' @export
remove.rows.any.nas <- function(x, indices=FALSE) {
    i <- plyr::aaply(x, 1, function(r) {any(is.na(r))})
    if (indices) {
        return(which(i))
    } else {
        return(x[ ! i, ])
    }
}

#' remove columns from a matrix that include any NAs
#'
#' @inheritParams remove.rows.all.nas
#' @return the matrix with columns removed, unless indices is TRUE
#' @export
remove.cols.all.nas <- function(x, indices=FALSE) {
    i <- plyr::aaply(x, 2, function(r) {all(is.na(r))})
    if (indices) {
        return(which(i))
    } else {
        return(x[ , ! i])
    }
}

#' remove columns from a matrix that are entirely NAs
#'
#' @inheritParams remove.rows.all.nas
#' @return the matrix with columns removed, unless indices is TRUE
#' @export
remove.cols.any.nas <- function(x, indices=FALSE) {
    i <- plyr::aaply(x, 2, function(r) {any(is.na(r))})
    if (indices) {
        return(which(i))
    } else {
        return(x[, ! i])
    }
}


add.missing.rows <- function(x, i) {
    i <- i[! i %in% rownames(x)]
    x <- rbind(x, matrix(NA, nrow=length(i), ncol=ncol(x),
                         dimnames=list(i, colnames(x))))
}

add.missing.cols <- function(x, j) {
    j <- j[! j %in% colnames(x)]
    x <- cbind(x, matrix(NA, nrow=nrow(x), ncol=length(j),
                         dimnames=list(rownames(x), j)))
}


#' align dataset margins
#'
#' @param datasets list of matrices (columns are samples, rows are genes)
#' @param na.rows if TRUE, add NAs for missing rows, otherwise they are removed
#' @param na.cols if TRUE, add NAs for missing columns, otherwise they are removed
#' @param use.dims either an integer or character vector of length 1
#' @return new list of matrices
#' @export
align.dataset.margins <- function(datasets, na.rows=TRUE, na.cols=FALSE,
                                  use.dims=NULL) {

    data.rownames <- plyr::llply(datasets, rownames)
    data.colnames <- plyr::llply(datasets, colnames)

    if (is.null(use.dims)) {
        if (na.rows) {
            row.collapse <- union
        } else {
            row.collapse <- intersect
        }

        if (na.cols) {
            col.collapse <- union
        } else {
            col.collapse <- intersect
        }

        data.rownames <- Reduce(row.collapse, data.rownames)
        data.colnames <- Reduce(col.collapse, data.colnames)

        if (length(data.rownames) < 1) stop("Error: no common rownames")
        if (length(data.colnames) < 1) stop("Error: no common colnames")

        if (na.rows) {
            datasets <- plyr::llply(datasets, add.missing.rows, data.rownames)
        }
        if (na.cols) {
            datasets <- plyr::llply(datasets, add.missing.cols, data.colnames)
        }
    } else {
        stopifnot(!is.null(data.rownames[[use.dims]]))
        stopifnot(!is.null(data.colnames[[use.dims]]))

        data.rownames <- data.rownames[[use.dims]]
        data.colnames <- data.colnames[[use.dims]]

        datasets <- plyr::llply(datasets, add.missing.rows, data.rownames)
        datasets <- plyr::llply(datasets, add.missing.cols, data.colnames)

    }

    datasets <- plyr::llply(datasets,
                      function(d) d[data.rownames, data.colnames])

}


#' make a tidy dataset
#' @inheritParams align.dataset.margins
#' @param dim.names character vector of length two -
#' these become the column labels in the tidy dataset
#' for the rownames and colnames of the original matrices
#' @return a tidy dataframe
#' @description Takes a list of matrices with similar row and column names
#' and runs tidyr::gather across each to make a "long" dataset
#' with a column of values for each matrix in the list.
#'
#' @export
make.a.tidy.dataset <- function(datasets, na.rows=TRUE, na.cols=FALSE,
                                use.dims=NULL,
                                dim.names=c("Gene", "Sample")) {

    datasets <- align.dataset.margins(datasets,
                                      na.rows=na.rows, na.cols=na.cols,
                                      use.dims=use.dims)

    tmp.df <- as.data.frame(datasets[[1]])

    rows <- rownames(tmp.df)
    cols <- colnames(tmp.df)

    tmp.df[, dim.names[1]] <- rows
    tmp.df <- tidyr::gather_(tmp.df, dim.names[2], "Temp", cols)
    tmp.df <- dplyr::select(tmp.df, -Temp)

    data.names <- names(datasets)
    names(data.names) <- data.names

    data.columns <- plyr::llply(data.names, function(d) {
        gath.df <- tidyr::gather_(as.data.frame(datasets[[d]]), dim.names[2], d, cols)
        return(gath.df[,d])
    })

    return(cbind(tmp.df, dplyr::as_data_frame(data.columns)))

}


#' convert data frame to matrix
#' @param df data frame
#' @return a matrix
#' @description Converts a data frame to a matrix
#' @export
df.to.mat <- function(df) {
    row.col <- colnames(df)[1]
    df %>%
        as.data.frame() %>%
        dplyr::select(1:3) %>%
        tidyr::spread_(colnames(df)[2], colnames(df)[3]) %>%
        set_rownames(.[[row.col]]) %>%
        dplyr::select_(paste("-", row.col)) %>%
        as.matrix()
}


#' convert matrix to data frame
#' @param mat matrix
#' @param row.name column label of resulting data frame from \code{mat} rows
#' @param col.name column label of resulting data frame from \code{mat} cols
#' @param dat.name column label of resulting data frame from \code{mat} data
#' @return a data frame
#' @description Converts a matrix to data frame
#' @export
mat.to.df <- function(mat, row.name="Row", col.name="Col", dat.name="Dat") {
    mat %>%
        as.data.frame %>%
        dplyr::mutate_(.dots=c(~rownames(.)) %>% set_names(row.name)) %>%
        tidyr::gather_(col.name, dat.name, colnames(mat)) %>%
        as.tbl
}


#' convert matrix to data frame
#' @param mat matrix
#' @param row.name column label of resulting data frame from \code{mat} rows
#' @param col.name column label of resulting data frame from \code{mat} cols
#' @param dat.name column label of resulting data frame from \code{mat} data
#' @return a data frame
#' @description Converts a matrix to data frame
#' @importFrom rlang enquo
#' @importFrom rlang "!!"
#' @importFrom plyr daply
#' @export
collapse_rows_of_matrix <- function(mat, group_df,
                                    collapse_fun=Matrix::colMeans,
                                    group_var=Group,
                                    sample_var=Sample,
                                    do_parallel=F, ...) {

    GroupVar <- enquo(group_var)
    SampleVar <- enquo(sample_var)

    GroupStr <- deparse(substitute(group_var))
    SampleStr <- deparse(substitute(sample_var))

    group_df <- group_df %>%
        dplyr::filter(!is.na(!!GroupVar), !is.na(!!SampleVar)) %>%
        dplyr::filter((!!SampleVar) %in% rownames(mat)) %>%
        dplyr::distinct(GroupVar, SampleVar)

    single_group_df <- group_df %>% dplyr::group_by(!!GroupVar) %>% dplyr::filter(n() == 1)
    group_df <- group_df %>% dplyr::group_by(!!GroupVar) %>% dplyr::filter(n() > 1)

    mat_1 <- mat[single_group_df[[SampleStr]],, drop=F] %>%
        set_rownames(single_group_df[[GroupStr]])

    mat_n <-
        daply(group_df, GroupVar, function(g) {
            m <- mat[unique(g[[SampleStr]]),, drop=F]
            do.call(collapse_fun, list(m, ...))
        }, .parallel=do_parallel)

    rbind(mat_1, mat_n) %>% {.[sort(rownames(.)), ]}

}

