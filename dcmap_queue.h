// DCMAP Queue class
#define QBLOCKSZ 1000

class dcmapq {
public:
	dcmapq() {}
	// Data members
	sp_mat SQ;								// sparse matrix for Q[entry i, branch/cluster/layer/ghat/G]
	sp_umat qinds;							// sparse vector of qinds (offset from first entry on SQ for that branch)
	uword b;								// pop: b
	uword k;								// pop: k
	uword l;								// pop: l
	double Ghat;							// pop: Ghat
	double G;								// pop: G
	uword Ncs;								// Ncs, use to index row of SQ given branch
	// Function members
	void initq(uword qsize, uword _Ncs);
	uword size();
	void addtoq(uword b, uword k, uword l, double Ghat, double G);
	void addtoq(uword b1, umat& qadd, double ghat1, double G1, double brpl);
	void popq(uword opt);
	uword peekq(uword opt);
	uword peekq(sp_imat& br_pl);
	void popg(uword opt, double Gmin);
	uword batchrm_bbrpl(uword branch, uword layer);
	uword batchrm_ggmin(double Gmin);
	uword batchmv_bbrpl(dcmapq& toQ, uword branch, uword layer);
	void batchcp_bbrpl(dcmapq& toQ, uword branch, uword layer, uword newbranch);
	mat exportq(string label);
};

// DCMAP Queue functions
// ----- initq
// Initialise sparse matrix Q with qsize = #rows, Ncs #rows per branch
void dcmapq::initq(uword Np, uword _Ncs) {
	SQ = sp_mat(Np*_Ncs, 5);
	Ncs = _Ncs;
	qinds = sp_umat(Np, 1);
}

// ----- size
// Returns size of Q
uword dcmapq::size() {
	return SQ.col(0).n_nonzero;
}

// ----- addtoq
// Add to Q, Qghat and Ql given branch, cluster, layer, Ghat and G
void dcmapq::addtoq(uword b1, uword k1, uword l1, double Ghat1, double G1) {
	rowvec q0 = { (double)b1, (double)k1,(double)l1, Ghat1, G1 };// rowvec to add
	SQ.row(Ncs * (b1-1) + qinds(b1 - 1, 0)) = q0;
	qinds(b1 - 1, 0)++;
	if (qinds(b1 - 1, 0) > Ncs) {
		cout << "=========================================================b,k,l,Ghat,G = "<< b1 << "," << k1 << "," << l1 << "," << Ghat1 << "," << G1 << "," << endl;
	}
}

// ----- exportq
// Export Q, Qghat and Ql (non-zero elements) to a dense matrix; pass label ---noprint to suppress printing.
mat dcmapq::exportq(string label) {
	//cout << (*this).size() << endl;
	mat M((*this).size(), 5);
	uword i = 0;
	for (auto it = SQ.col(0).begin(); it != SQ.col(0).end(); ++it) {
		for (uword j = 0; j < 5; j++) {
			//cout << it.row() << "," << i << endl;
			M(i, j) = SQ(it.row(), j);
		}
		i++;
		//cout << endl;
	}
	if (label.compare("---noprint") != 0)				// compare returns 0 if strings are equal
		M.print(label);
	return M;
}

// ----- popq
// Pop Q using key opt, 0=ghat, 1=l. Save results to members b,k,l,Ghat,G. 
void dcmapq::popq(uword opt) {
	uword i = (*this).peekq(opt);
	SQ.row(i).clean(arma::math::inf());					// set to zero the popped row in Q sparse mat
}

// ----- peekq
// Peek Q using key opt, 0=ghat, 1=l. Save results to members b,k,l,Ghat,G. 
// Returns i = row of the min value
uword dcmapq::peekq(uword opt) {
	uword i=0, col;										// i to store min row, col of SQ
	double val = arma::math::inf();						// to store min Ghat or l
	if (opt == 0) {										// use ghat min, col3
		col = 3;
	}
	else {
		col = 2;										// use l min, col2
	}

	auto it = min_element(SQ.col(col).begin(), SQ.col(col).end());
	val = *it;
	i = it.row();

	//for (auto it = SQ.col(0).begin(); it != SQ.col(0).end(); ++it) {// foreach non-zero val
	//	if (SQ(it.row(), col) < val) {
	//		val = SQ(it.row(), col);
	//		i = it.row();
	//	}
	//}
	b = (uword)SQ(i, 0);								// write popped value to b,k,l,Ghat,G
	k = (uword)SQ(i, 1);
	l = (uword)SQ(i, 2);
	Ghat = SQ(i, 3);
	G = SQ(i, 4);
	return i;
}


// ----- peekq
// Peek Q using key opt, 0=ghat, 1=l. Save results to members b,k,l,Ghat,G. Looks for first
// Returns i = row of the min value
uword dcmapq::peekq(sp_imat& br_pl) {
	uword i = 0;										// i to store min row
	double val = arma::math::inf();						// to store min Ghat or l
	
	auto itb = SQ.col(0).begin();
	auto ite = SQ.col(0).end();

	//cout << (--ite).row() << " vs " << next(itb, SQ.col(0).n_nonzero-1).row() << endl;

	//if (SQ.col(0).n_nonzero > 2 * QBLOCKSZ) {
	//	double dix = randu<double>(distr_param(0, SQ.col(0).n_nonzero));	// start index as offset from SQ.col(col).begin()
	//	uword diy = min((uword)dix + 2 * QBLOCKSZ, SQ.col(0).n_nonzero);
	//	//cout << dix << "---" << diy << endl;
	//	itb = next(itb, (uword)max(dix - QBLOCKSZ, 0.0));
	//	ite = next(SQ.col(0).begin(), diy);
	//	//cout << itb.row() << "," << ite.row() << " out of " << SQ.col(0).n_nonzero << endl;
	//}

	for (auto it = itb; it != ite; ++it) {// foreach non-zero val
		if (SQ(it.row(), 3) < val && SQ(it.row(), 2) <= (double)abs(br_pl(SQ(it.row(),0)-1, 0))) {
			val = SQ(it.row(), 3);
			i = it.row();
		}
	}
	b = (uword)SQ(i, 0);								// write popped value to b,k,l,Ghat,G
	k = (uword)SQ(i, 1);
	l = (uword)SQ(i, 2);
	Ghat = SQ(i, 3);
	G = SQ(i, 4);
	return i;
}


// ----- popq
// Pop Q using key opt, 0=ghat, 1=l. Save results to members b,k,l,Ghat,G. Keep popping until G<=Gmin
void dcmapq::popg(uword opt, double Gmin) {
	G = datum::inf;
	while (G > Gmin && (*this).size()>0) {
		(*this).popq(opt);
	}
}


// ----- batchrm_bbrpl
// Batch remove from *this dcmap queue entries matching branch b and l <= layer (=brpl(b)); returns #removed
// Set layer >= lmax+1 to remove all entries on a branch
uword dcmapq::batchrm_bbrpl(uword branch, uword layer) {
	uword nummved = 0;									// number removed from *this queue
	uword sqind = Ncs * (branch - 1);					// ind to start of SQ rows for branch
	for (uword i = sqind; i < sqind + qinds(branch - 1,0); i++) {// for Q entries up to qinds in that branch
		if (SQ(i, 0) != 0 && SQ(i, 2) <= (double)layer) {
			SQ.row(i).clean(arma::math::inf());
			nummved++;
		}
	}
	return nummved;
}

// ----- batchmv_bbrpl
// Batch move from *this dcmap queue entries matching branch and l <= layer (=brpl(b)) to another dcmap queue
uword dcmapq::batchmv_bbrpl(dcmapq& toQ, uword branch, uword layer) {
	uword nummved = 0;									// number removed from *this queue
	uword sqind = Ncs * (branch - 1);					// ind to start of SQ rows for branch
	for (uword i = sqind; i < sqind + qinds(branch - 1, 0); i++) {// for Q entries up to qinds in that branch
		if (SQ(i, 0) != 0 && SQ(i, 2) <= (double)layer) {		// if Q(layer) <= brpl(b)
			toQ.addtoq(branch,SQ(i,1),(uword)SQ(i,2),SQ(i,3),SQ(i,4));
			nummved++;
			SQ.row(i).clean(arma::math::inf());			// also delete from source Q
		}
	}
	return nummved;
}


// ----- batchcp_bbrpl
// Batch copy from *this dcmap queue entries matching branch and l == layer (=brpl(b)) to another dcmap queue
void dcmapq::batchcp_bbrpl(dcmapq& toQ, uword branch, uword layer, uword newbranch) {
	//uword nummved = 0;									// number copied from *this queue
	uword sqind = Ncs * (branch - 1);					// ind to start of SQ rows for branch
	for (uword i = sqind; i < sqind + qinds(branch - 1, 0); i++) {// for Q entries up to qinds in that branch
		if (SQ(i, 2) == (double)layer && SQ(i, 0) != 0) {		// if Q(layer) == brpl(b)
			toQ.addtoq(newbranch, SQ(i, 1), (uword)SQ(i, 2), SQ(i, 3), SQ(i, 4));
			//nummved++;
		}
	}
	//return nummved;
}


// ----- addtoq
// Add to the queue entries on branch b1, k-l in vector<vector<double> >qadd (nx2 matrix) 
// ghat1 and G1 values, brpl is brpl(b) - for Qdash, this checks each qadd row to addtoq; for Q, set brpl=lmax+1 so always true
// if it exists on Q, update Q if proposed ghat is lower; else add to Q
void dcmapq::addtoq(uword b1, umat& qadd, double ghat1, double G1, double brpl) {
	uword sqind = Ncs * (b1 - 1);						// ind to start of SQ rows for branch
	// Update existing entries
	for (uword i = sqind; i < sqind + qinds(b1 - 1, 0); i++) {// for Q entries up to qinds in that branch
		if (SQ(i, 0) != 0.0 && SQ(i, 2) <= brpl) {		// if Q(layer) <= brpl(b)
			for (uword j = 0; j < qadd.n_rows; j++) {	// foreach entry to add / update
				if (qadd.row(j)(0) == (uword)SQ(i, 1) && qadd.row(j)(1) == (uword)SQ(i, 2)) {// match for k-l (and b before)
					if (ghat1 < SQ(i, 3)) {				// proposed entry ghat < existing entry - update existing
						SQ(i, 3) = ghat1;
						SQ(i, 4) = G1;
					} // end if update existing entry
					qadd.row(j)(0) = 0;					// mark as already updated in Q
				}
			} // end foreach entry to add
		} // end if active branch brpl
	} // end foreach existing branch entry
	
	for (uword i = 0; i < qadd.n_rows; i++) {// foreach entry to add 
		if (qadd.row(i)(0) > 0 && (double)qadd.row(i)(1) <= brpl) {	// not already updated in Q *and* eligible to add/pop in the case of Qdash
			(*this).addtoq(b1, qadd.row(i)(0), qadd.row(i)(1), ghat1, G1);
			qadd.row(i)(0) = 0;							// mark as already updated in Q
		}
	}
}


// ----- batchrm_ggmin
// Batch remove from *this dcmap queue entries where G > Gmin; returns #removed
uword dcmapq::batchrm_ggmin(double Gmin) {
	uword nummved = 0;									// number removed from *this queue
	for (auto it = SQ.col(0).begin(); it != SQ.col(0).end(); ++it) {// foreach non-zero val
		if (SQ(it.row(), 4) > Gmin) {
			SQ.row(it.row()).clean(arma::math::inf());
		}
	}
	return nummved;
}



//// find current non-zero entries in SQ on this branch
//auto it0 = SQ.col(0).begin();
//auto it1 = SQ.col(0).end();
//uword nnz = SQ.col(0).n_nonzero;						// number of non-zero SQ entries
//uword ilow = (b1 - 1) * Ncs;							// >= lower bound SQ entries for branch b1
//uword iupp = b1 * Ncs;								// < upper bound SQ entries for branch b1




