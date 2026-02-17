getMM_new <- function(object, vc = NULL, recordsToUse = NULL) {
  if (is.null(vc)) {
    vc <- VarCorr(object)
  }
  n <- length(object@resp$y)
  vc_e <- attr(VarCorr(object), "sc")^2
  Ri <- Matrix::Diagonal(n) * (1/vc_e)
  X <- getME(object, "X")
  Z <- getME(object, "Z")
  y <- getME(object, "y")
  if (!is.null(recordsToUse)) {
    X <- X[recordsToUse, , drop = FALSE]
    Z <- Z[recordsToUse, , drop = FALSE]
    y <- y[recordsToUse]
    Ri <- Ri[recordsToUse, recordsToUse]
  }
  C11 <- t(X) %*% Ri %*% X
  C12 <- t(X) %*% Ri %*% Z
  C21 <- t(Z) %*% Ri %*% X
  C22 <- t(Z) %*% Ri %*% Z
  C <- rbind(cbind(C11, C12), cbind(C21, C22))
  nm <- c(colnames(X), colnames(Z))
  colnames(C) <- rownames(C) <- nm
  RHS1 <- t(X) %*% Ri %*% y
  RHS2 <- t(Z) %*% Ri %*% y
  RHS <- rbind(RHS1, RHS2)
  fl <- object@flist
  asgn <- attr(fl, "assign")
  pnms <- names(object@flist)
  Gi <- Matrix::Matrix(0, nrow = nrow(C), ncol = ncol(C))
  for (iFac in pnms) {
    tn <- which(match(iFac, names(fl)) == asgn)
    vcov <- do.call(Matrix::bdiag, vc[tn])
    vcov <- vcov + diag(1e-06, ncol(vcov), ncol(vcov))
    LLt <- Matrix::Diagonal(length(unique(object@flist[[iFac]])))
    rowsi <- list()
    for (j in 1:length(tn)) {
      ind <- (object@Gp)[tn[j]:(tn[j] + 1L)]
      # rowsi[[j]] <- ((ind[1] + 1L):ind[2]) + 1  # Error
      rowsi[[j]] <- ((ind[1] + 1L):ind[2]) + ncol(X)
    }
    Gi[unlist(rowsi), unlist(rowsi)] <- kronecker(LLt, solve(Matrix::nearPD(vcov)$mat))
  }
  C <- C + Gi + diag(1e-04, ncol(C), ncol(C))
  C_inv <- solve(C)
  rownames(C_inv) <- colnames(C_inv) <- c(colnames(X), colnames(Z))
  bu <- C_inv %*% RHS
  rownames(bu) <- rownames(C_inv)
  if (length(object@relfac) > 0) {
    ROT <- Matrix::Diagonal(n = nrow(C))
    for (iFac in pnms) {
      tn <- which(match(iFac, names(fl)) == asgn)
      for (j in 1:length(tn)) {
        ind <- (object@Gp)[tn[j]:(tn[j] + 1L)]
        rowsi <- ((ind[1] + 1L):ind[2]) + ncol(X)
        if (iFac %in% names(object@relfac)) {
          pick <- rownames(C)[rowsi]
          ROT[rowsi, rowsi] <- object@relfac[[iFac]][pick,
                                                     pick]
        }
      }
    }
    rn <- rownames(C_inv)
    buROT <- t(as.matrix(t(bu) %*% ROT))
    C_invROT <- t(ROT) %*% C_inv %*% (ROT)
    rownames(buROT) <- rn
    colnames(C_invROT) <- rownames(C_invROT) <- rn
    return(list(Ci = C_invROT, bu = buROT, RHS = RHS))
  }
  else {
    return(list(Ci = C_inv, bu = bu, RHS = RHS))
  }
}
