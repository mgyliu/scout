# Location and scale comparison: default vs. cellwise

X_mat <- matrix(rnorm(n=100), ncol=4) 
Y_mat <- as.matrix(rnorm(n=25), ncol=1)

mean(Y_mat); mean_cw(Y_mat)
sd(Y_mat); scale_cw(Y_mat)

apply(X_mat, 2, mean)
apply(X_mat, 2, mean_cw)

apply(X_mat, 2, sd)
apply(X_mat, 2, scale_cw)

scaledX_1 <- scale(X_mat, center=TRUE, scale=TRUE)
scaledX_2 <- scale(X_mat, center=apply(X_mat, 2, mean), scale=apply(X_mat, 2, sd))

all(scaledX_1 == scaledX_2)

scaledX_3 <- scale(X_mat, center=apply(X_mat, 2, mean_cw), scale=apply(X_mat, 2, scale_cw))

scaledX_3 - scaledX_1
