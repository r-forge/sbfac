.givensrot <- function(t, dim=2, i=1, j=2) {
	rot = matrix(rep(0, dim^2), nrow=dim)
	rot[i,i] = cos(t)
	rot[j,j] = cos(t)
	rot[i,j] = sin(t)
	rot[j,i] = -sin(t)
	return(rot)
}