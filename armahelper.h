#pragma once

// Global functions
//#define DOUBLEPRECISION 0.0001
uvec mod(uvec a, uword n) {										// modulus operator for uvec modulo scalar n
	return a - arma::floor(a / n) * n;
}

uvec std_setdiff(const uvec& x, const uvec& y) {				// setdiff of two uvecs
	// from https://stackoverflow.com/questions/29724083/trying-to-write-a-setdiff-function-using-rcpparmadillo-gives-compilation-error
	std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
	std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
	std::vector<int> out;

	std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
		std::inserter(out, out.end()));

	return arma::conv_to<arma::uvec>::from(out);
}

template <typename t>
std::vector<std::vector<t> > GetUniqueRows(std::vector<std::vector<t> > input) {// get unique rows vector<vector<>>
	// from https://stackoverflow.com/questions/3169960/determining-the-unique-rows-of-a-2d-array-vectorvectort
	std::sort(input.begin(), input.end());
	input.erase(std::unique(input.begin(), input.end()), input.end());
	return input;
}

template <typename t>
Mat<t> GetNondupRows(Mat<t>& x) { // return matrix non-duplicated rows of matrix x
	// code from: https://stackoverflow.com/questions/37143283/finding-unique-rows-in-armamat
	uvec mask(x.n_rows, fill::zeros);
	for (uword i = 0; i < x.n_rows; i++) {
		for (uword j = i + 1; j < x.n_rows; j++) {
			if (all(x.row(i) == x.row(j))) {
				mask(j) = 1;
				break;
			}
			//if (approx_equal(conv_to<rowvec>::from(x.row(i)), conv_to<rowvec>::from(x.row(j)), "absdiff", DOUBLEPRECISION)) {
			//	mask(j) = 1;
			//	break;
			//}
		}
	}
	return x.rows(find(mask == 0));
}

template <typename t>
Col<t> GetNondupVec(Col<t>& x) { // return vector non-duplicated elements of column vector x
	// code from: https://stackoverflow.com/questions/37143283/finding-unique-rows-in-armamat
	uvec mask(x.size(), fill::zeros);
	for (uword i = 0; i < x.size(); i++) {
		for (uword j = i + 1; j < x.size(); j++) {
			if (x(i) == x(j)) {
				mask(j) = 1;
				break;
			}
		}
	}
	return x(find(mask == 0));
}

template <typename t>
Col<t> GetNondupCol(Col<t>& x) { // return column vector non-duplicated entries of vector x
	// code from: https://stackoverflow.com/questions/37143283/finding-unique-rows-in-armamat
	uvec mask(x.size(), fill::zeros);
	for (uword i = 0; i < x.size(); i++) {
		for (uword j = i + 1; j < x.size(); j++) {
			if (x(i) == x(j)) {
				mask(j) = 1;
				break;
			}
		}
	}
	return x(find(mask == 0));
}

//template <typename t>
//Mat<t> GetNondupRows(Mat<t>& x) { // return matrix non-duplicated rows of matrix x
//	// code from: https://stackoverflow.com/questions/37143283/finding-unique-rows-in-armamat
//	uword count = 1, nr = x.n_rows, nc = x.n_cols;
//	Mat<t> result(nr, nc);
//	result.row(0) = x.row(0);
//
//	for (uword i=1; i < nr; i++) {
//		bool matched = false;
//		if (approx_equal(x.row(i), result.row(0), "absdiff", DOUBLEPRECISION)) continue;
//
//		for (uword j = i + 1; j < nr; j++) {
//			if (approx_equal(x.row(i), x.row(j),"absdiff",DOUBLEPRECISION)) {
//				matched = true;
//				break;
//			}
//		}
//		if (!matched) result.row(count++) = x.row(i);
//	}
//	return result.rows(0, count - 1);
//}

template <typename t>
void setspmscal(SpMat<t>& X, const uvec& rowinds, const uvec& colinds, t scal) {	// set X(rowinds,colinds) = scalar
	for (uword i = 0; i < rowinds.size(); i++) {
		for (uword j = 0; j < colinds.size(); j++) {
			X(rowinds(i), colinds(j)) = scal;
		}
	}
}

template <typename t>
void setspmscal(SpMat<t>& X, const uvec& rowinds, const uword colind, t scal) {	// set X(rowinds,col) = scalar
	for (uword i = 0; i < rowinds.size(); i++) {
		X(rowinds(i), colind) = scal;
	}
}

template <typename t>
void setspmvec(SpMat<t>& X, const uvec& rowinds, const uword colind, const Col<t>& Y) {	// set X(rowinds,col) = vector
	for (uword i = 0; i < rowinds.size(); i++) {
		X(rowinds(i), colind) = Y(i);
	}
}

template <typename t>
vector<uword> findspcol(SpSubview_col<t>& X, t y) {						// find uword inds to sparse matrix column vec X equal to y
	vector<uword> inds;
	for (auto it = X.begin(); it != X.end(); ++it) {
		if (*it == y)
			inds.push_back(it.row());
	}
	return inds;
}

//template <typename t>
//vector<uword> findspcol_gt(SpSubview_col<t>& X, t y) {					// find uword inds to sparse matrix column vec X > y
//	vector<uword> inds;
//	for (auto it = X.begin(); it != X.end(); ++it) {
//		cout << *it << ", " << y << endl;
//		if (*it > y)
//			inds.push_back(it.row());
//	}
//	return inds;
//}

template <typename t>
vector<uword> findsprow(SpSubview_row<t>& X, t y) {						// find uword inds to sparse matrix row vec X equal to y
	vector<uword> inds;
	for (auto it = X.begin(); it != X.end(); ++it) {
		if (*it == y)
			inds.push_back(it.col());
	}
	return inds;
}

template <typename t>
bool findspcol_anyneq(SpSubview_col<t>& X, uvec& subinds, t y) {		// return TRUE/FALSE if any in sparse matrix column vec X with subinds *not* equal to y
	for (auto it = X.begin(); it != X.end(); ++it) {
		if (any(subinds == it.row()) && *it != y)
			return true;
	}
	return false;
}

template <typename t>
bool findspcol_allnz(SpSubview_col<t>& X, uvec& subinds) {				// return TRUE/FALSE if all in sparse matrix column vec X with subinds are non-zero
	uword i = 0;
	for (auto it = X.begin(); it != X.end(); ++it) {
		if (any(subinds == it.row()))
			i++;
	}
	if (i == subinds.size()) {
		return true;
	}
	else {
		return false;
	}
}

// find uword inds to sparse matrix row vec X that satisfies condition type (0=equal,1=greater) to scalar y for indices of X subinds
template <typename t>
vector<uword> findsprow(SpSubview_row<t>& X, t y, uvec& subinds, uword type) {
	vector<uword> inds;
	for (auto it = X.begin(); it != X.end(); ++it) {
		if (type == 1 && any(subinds == it.col()) && *it > y) {
			inds.push_back(it.col());
		}
		else if (type == 0 && any(subinds == it.col()) && *it == y) {
			inds.push_back(it.col());
		}
	}
	return inds;
}

// find min non-zero element in sparse matrix row vec X for indices of X subinds
template <typename t>
t minsprow(SpSubview_row<t>& X, uvec& subinds) {
	uword ind;
	t val = datum::inf;
	for (auto it = X.begin(); it != X.end(); ++it) {
		if (any(subinds == it.col()) && (*it) < val) {
			val = *it;
			ind = it.col();
		}
	}
	return val;
}

template <typename t>
vector<t> spcolsubset(SpSubview_col<t>& X, uvec& subinds) {						// return subset of X as full vector given uvec indices in inds
	vector<t> Y;
	for (auto it = X.begin(); it != X.end(); ++it) {
		if (any(subinds == it.row()))
			Y.push_back(*it);
	}
	return Y;
}
//// find uword inds to sparse matrix 2 columns: vec X1 equal to y1, X2 ==y2; assumes X1.size()==X2.size() and same type
//template <typename t>
//vector<uword> findspcol(SpSubview_col<t>& X1, SpSubview_col<t>& X2, t y1, t y2) {
//	vector<uword> inds;
//	for (auto it = X1.begin(); it != X1.end(); ++it) {
//		if (*it == y1 && X2[it.row()] == y2)
//			inds.push_back(it.row());
//	}
//	return inds;
//}