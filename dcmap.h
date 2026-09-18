

// Includes
#include <armahelper.h>
#include <dcmap_cost.h>
#include <dcmap_queue.h>


// DCMAP class
class dcmap {
public:
	dcmap(uword& _numnodes, uvec& _numstates, uvec& _leaves, uvec& _roots, mat& _dag, uword _maxiter);
	// BN Data members
	uword numnodes;						// number of nodes
	uvec numstates;						// uvec[x]=number of states
	uvec leaves;						// vec of leaf node indices
	uvec roots;							// vec of root node indices
	mat dag;							// DAG[i,j] = datum::nan (no arc), t>=0 for time slice t of parent j to node i
	
	// Data members
	uvec cl_layer;						// layer foreach node
	uvec cl_order;						// cl_order[x]=order, x=node
	uword lmax;							// max layer
	uword prodps, Np, Ncs;				// max branching size Np; prodps = search space size; sparse queue size per branch
	sp_umat Uhat;						// Uhat[k-i,b] = {0,1}
	sp_umat U;							// U[i,b] = k
	sp_mat J;							// J[k-l-i,b] = joint partial result
	sp_mat G;							// G[l,b] = objective function
	double Gmin;						// Gmin
	sp_umat br_lk;						// branch linking via br_lk and br_pl: br_lk[b_il,b_jl]=1 if linked at layer l
	sp_imat br_pl;						// branch popped layer (for queue peek)
	dcmapq Q;							// DCMAP Queue object
	dcmapq Qdash;						// DCMAP Queue object for entries eligible to be popped br_pl(b)>=q_b(l)
	uword bmax = 1;						// next branch index to use
	umat sameclust;						// sameclust matrix structure: col=x, row=nodes that can be in same clust as x, 1=yes/0=no

	uword iter;							// number of iterations
	uword maxiter;						// max # iterations
	mat debugdat;						// debugging matrix[iter-b-k-l-G.l-ghat-Gmin-delb-node0-node1-..., iteration]
	uword debugi;						// index to debugdat
	uword Gminind;						// last iteration when Gmin was updated

	// Function members	
	void get_cl_layer();
	void get_branching_size();
	void get_cl_order();
	void makesameclust();
	uvec get_intlk(const uvec& xind, const uword bdash, const uword k);				// returns 1=int node, 2=lk node
	void search(double seed);			// DCMAP search
	void search_breadth(double seed);	// breadth first search for comparison
	// Note that search() depends on cost functions defined in dcmap_cost.h of: propchild, bwdselimcost, heuristic
};

// DCMAP functions
// ----- get_cl_layer
// Populates cl_layer and lmax data members given dag, numnodes and leaves
void dcmap::get_cl_layer() {
	uword writeind = 0;										// layer index to write
	uvec inds(leaves);										// node inds to write to
	cl_layer.resize(numnodes);								// initialise cl_layer

	while (inds.size() > 0) {								// while there are unassinged nodes
		uvec writevec(inds.size(), fill::value(writeind));	// turn writeind into a vector
		cl_layer(inds) = max(writevec, cl_layer(inds));		// write max
		writeind++;
		// get parents of nodes in inds to update next layer (find returns columnwise inds, so need modulo)
		inds = unique( mod(find((dag.rows(inds) >= 0).t()),numnodes) );
	}
	lmax = writeind-1;
}


// ----- get_branching_size
// Populates Np as an estimate of max branching size given cl_layer, dag
void dcmap::get_branching_size() {
	uvec ps(lmax + 1, fill::value(1));						// layer 0 maps to ps[0]; ps[0] = 1
	uvec ps2(ps);											// prod of consecutive layers of ps
	uword ncs = sum((cl_layer == 0)); 						// cumsum for sparse queue size foreach branch	- for layer i
	Ncs = ncs;												//												- total cumsum
	for (uword i = 1; i < ps.size(); i++) {
		ps(i) = prod(sum(dag.cols(find(cl_layer == i))>=0,0) + 1);
		ps2(i) = ps(i) * ps(i - 1);
		ncs = sum((cl_layer == i)) + ncs;
		Ncs += ncs;
	}
	Np = ps2.max();											// max branching size
	prodps = prod(ps);										// search space size
	if (Np < 10000) {
		Np *= 4;
	} else if(Np > 100000) {
		Np /= 20;
	}
	//Ncs *= 10;
	//Np *= 10;
}


// ----- get_cl_order
// Order of node elimination with layer by layer backwards elimination
void dcmap::get_cl_order() {
	cl_order = cl_layer;									// cl_order[x]=order, x=node
	// Use backwards order by layer from leaf nodes; unaffected by cluster def
	uword writeind = 0;										// order index to write
	vec cpt_cost(cl_layer.size());							// use to differentiate nodes within a layer
	for (uword i = 0; i < cl_layer.size(); i++) {
		cpt_cost(i) = numstates(i) * prod(numstates(find(dag.row(i) >= 0)));
	}
	for (uword i = 0; i <= max(cl_layer); i++) {
		uvec inds(find(cl_layer == i));
		cl_order(inds) = writeind + sort_index(cpt_cost(inds));
		writeind += inds.size();
	}
}

// ----- init
// Initialises dcmap object given numnodes, vector numstates foreach node, vec leaves (indexing numstates), vec roots (ditto)
// and matrix dag, mi=max #iterations
dcmap::dcmap(uword& _numnodes, uvec& _numstates, uvec& _leaves, uvec& _roots, mat& _dag, uword _maxiter) {
	// Copy BN data
	numnodes = _numnodes;
	numstates = _numstates;
	leaves = _leaves;
	roots = _roots;
	dag = _dag;

	get_cl_layer();											// init cl_layer
	get_branching_size();									// get Np and prodps
	get_cl_order();											// get cl_order
	leaves = leaves(sort_index(cl_order(leaves)));			// sort leaves

	Q.initq(Np,Ncs);										// init queue Q object
	Qdash.initq(Np,Ncs);									// init queue Qdash (eligible to be popped) object
	Uhat = sp_umat(numnodes*numnodes, Np);					// Uhat[k-i,b]
	U = sp_umat(numnodes, Np);								// U[i,b]
	J = sp_mat(numnodes * numnodes * (lmax+1), Np);				// J[k-l-i,b]
	G = sp_mat((lmax + 1), Np);								// G[l,b]
	br_lk = sp_umat(Np * (lmax + 1), Np * (lmax + 1));		// il=branch-layer, jl=branch-layer, y_ij=1 branches linked at layer l (allow links at multiple layers)
	br_pl = sp_imat(Np, 1);									// branch popped layer (for queue peek)
	bmax = 1;
	
	makesameclust();										// make sameclust matrix

	maxiter = _maxiter;
	debugdat = mat(8+numnodes,maxiter*4);					// debugging matrix[iter-l-k-b-G.l-ghat-Gmin-delb-node0-node1-..., iteration]
	debugi = 0;												// index to debugdat
}

// ----- makesameclust
// Creates sameclust matrix give numnodes and dag matrix
// col=x, row=nodes that can be in same clust as x, 1=yes/0=no
// row: given a node, what other nodes are connected by a parind
// col: when proposing clusts for parind, which nodes to propose similar clusts to
void dcmap::makesameclust() {
	sameclust = umat(numnodes, numnodes);
	vector<uvec> parinds;									// indices to parent nodes foreach node
	vector<uvec> pardesc(numnodes);							// parind + descendants down to leaf nodes

	for (uword i = 0; i < numnodes; i++) {
		parinds.push_back( unique(mod(find((dag.row(i) >= 0).t()), numnodes)) );
		// parinds.at(i).print("- ");
	}
	for (uword i = 0; i < numnodes; i++) {					// get all descendants of parinds of x and filter by layer ==x
		uvec inds = parinds.at(i);							// start at parent inds
		while (inds.size() > 0) {
			pardesc.at(i) = join_cols(pardesc.at(i), inds);	// add to pardesc
			inds = unique(mod(find((dag.cols(inds) >= 0)), numnodes));	// get children 
		} // end while
		//if (i == 12) {
		//	parinds.at(i).print("parinds ");
		//	pardesc.at(i).print("pardesc ");
		//}
		pardesc.at(i) = pardesc.at(i)(find(cl_layer(pardesc.at(i)) == cl_layer(i)));	//filter descendants by layer ==x
		pardesc.at(i) = std_setdiff(unique(pardesc.at(i)), uvec(1, fill::value(i)));
		//if(i==12)
		//	pardesc.at(i).print("pardesc ");
		uvec shedi;
		for (auto it = parinds.at(i).begin(); it != parinds.at(i).end(); ++it) {
			//cout << "it " << *it << endl;
			if (cl_layer((*it)) - cl_layer(i) > 1) {		// l(pardesc) > l(i)+1; check no path with length > 2
				inds = uvec(1,fill::value(*it));
				uvec pathlen(numnodes, fill::zeros);
				uvec dddd;
				uword xpathlen = 0;							// current path length
				while (inds.size() > 0) {
					inds = unique(mod(find((dag.cols(inds) >= 0)), numnodes));
					xpathlen++;
					for (auto jt = inds.begin(); jt != inds.end(); ++jt) {
						if (pathlen(*jt) < xpathlen)
							pathlen(*jt) = xpathlen;
					}
				} // end while
				if (pathlen(i) > 1) {						// if path length at node i > 1, orignating from *it=pardesc.at(i), rm it
					//pardesc.at(i).shed_row(*it);
					shedi = join_cols(shedi, uvec(1, fill::value(it-parinds.at(i).begin())));
				}
			}
		} // end filter out pardesc.at(i)
		parinds.at(i).shed_rows(shedi);

		urowvec x = sameclust.row(i);
		x(pardesc.at(i)).fill(1);
		x(parinds.at(i)).fill(1);
		sameclust.row(i) = x;
	} // end for
	//sameclust.print("sameclust ");
}

// ----- intlk
// returns 1=int node, 2=lk node; uses dcmap::dag, and dcmap::U with inputs xind uvec and bdash uword (current branch), current cluster k
uvec dcmap::get_intlk(const uvec& xind, const uword bdash, const uword k) {
	uvec intlk(xind.size(), fill::value(1));
	//cout << "k = " << k << ", b = " << bdash-1 << endl;
	//U.print("U: ");
	for (auto it = xind.begin(); it != xind.end(); ++it) {
		uvec inds = find((dag.col(*it) >= 0));	// get children 
		//inds.print("children: ");
		if (inds.size()>0 && findspcol_anyneq(U.col(bdash - 1), inds, k))
				intlk(it - xind.begin()) = 2;
	}
	return intlk;
}


// ----- search
// Run DCMAP search algorithm 
// Inputs: seed is the random seed to use
// Note that this depends on cost functions defined in dcmap_cost.h of: propchild, bwdselimcost, heuristic
void dcmap::search(double seed) {
	arma_rng::set_seed(seed);

	SpSubview_col<double> x = J.col(0);
	Gmin = heuristic(linspace<uvec>(0, numnodes - 1, numnodes), dag, uvec(U.col(0)), numnodes, 
		cl_layer, cl_order, numstates, x);											// get heuristic and update Gmin
	Gminind = 1;
	leaves = leaves(sort_index(cl_order(leaves)));									// order leaf nodes; set Uhat(leaf nodes, branch1) to 1
	setspmscal(Uhat, leaves + linspace<uvec>(0, leaves.size() - 1, leaves.size()) * numnodes, uvec(1,fill::zeros), (uword)1);
	for (uword i = 0; i < leaves.size(); i++) 										// add leaf nodes to Q
		Qdash.addtoq(1, i+1, 0, Gmin, 0);

	iter = 1;
	//uword prevpeek = 1;															// previous peek: 1 = popped Qdash or lowest layer entry, 0 = least cost entry
	while (iter <= maxiter){														// main loop
		if(iter % 50 ==0)
			cout << "Iteration " << iter << " Q size " << Q.size() << " Gmin " << Gmin << endl;
		// ----- Pop Qdash with random switching between min Ghat and min layer
		if (Qdash.size() == 0 && Q.size()>0) {										// Qdash empty (no active branches)
			//cout << "Gmin: " << Gmin << endl;
			//Q.batchrm_ggmin(Gmin);												// flush any dominated Q entries
			//uword delentry = 1;
			if (randu() > 0.6){//} || prevpeek == 0) {								// Pick lowest layer branch on Q to activate br.pl
				Q.peekq(1);
				//prevpeek = 1;
			}
			else {																	// Pick least cost branch b on Q to activate branch (br.pl)
				//Q.peekq(0);
				//prevpeek = 0;
				//delentry = 0;
				if (Q.peekq(br_pl) == 0)
					break;
			}
			//cout << "*" << Q.b << "," << Q.k << "," << Q.l << "-" << abs(br_pl(Q.b - 1)) << "," << Q.Ghat << "," << Q.G << endl;
			br_pl(Q.b - 1) = abs(br_pl(Q.b - 1));									// activate using br_pl and add to Qdash
			Q.batchmv_bbrpl(Qdash, Q.b, (uword)br_pl(Q.b - 1));						// batch move all activated entries Q to Qdash

		} else if(Qdash.size()>0) {													// Qdash not empty
			//prevpeek = 1;
			Qdash.popg(0,Gmin);														// Pop until G<Gmin
			//cout << Qdash.G << ", " << Gmin << endl;
			//cout << Qdash.b  << "," << Qdash.k << "," << Qdash.l << "," << Qdash.Ghat << ","  << Qdash.G  << endl;
			//cout << "-------------------- Iter " << iter << " b " << Qdash.b << " k " << Qdash.k << " l " << Qdash.l << endl;
			// ----- Prune dominated branches
			bool bdashprune = false;
			if (Qdash.l > 0) {														// pop layer l if U_{ l - 1 } not dominated given U^* _{ l - 2 }
				uvec bvec2(findsprow(br_lk.row(Qdash.l - 1 + (Qdash.b - 1) * (lmax + 1)), (uword)1));// linked branches layer l-1, adjust for index start at 0
				uvec bvec(bvec2.size() + 1);				
				bvec(0) = Qdash.b-1;	
				// convert bvec2 ind to br_lk to bvec ind to G
				if (bvec2.size() > 0) {												// only prune if there are linked branches
					bvec(span(1, bvec.size() - 1)) = (bvec2 + 1 - Qdash.l) / (lmax + 1);// adjust for index start at 0
					double Gbmin = minsprow(G.row(Qdash.l - 1), bvec);				// min G at layer l-1
					uvec inds(findsprow(G.row(Qdash.l-1), min(Gmin, Gbmin), bvec, 1));// dominated by Gmin || dominated by another branch at layer l-1

					for (auto it = inds.begin(); it != inds.end(); ++it) {			// prune - del from br_pl and Q, Qdash
						br_pl(*it) = 0;
						Q.batchrm_bbrpl(*it + 1, lmax + 1);							// br_pl index starts at 0, so branch is ind+1
						Qdash.batchrm_bbrpl(*it + 1, lmax + 1);
					}
					debugdat(7, debugi) = inds.size();								// debugdat - store number pruned branches
					if (any(inds == (Qdash.b - 1))) {								// popped branch was pruned; indexing starts at 0, hence b-1
						bdashprune = true;
						debugdat(span(0, 3), debugi) = { (double)iter, (double)Qdash.b, (double)Qdash.k, (double)Qdash.l };// update debugdat
						debugi++;
					}
				} // end if bvec2 not empty
			} // end prune dominated branches
			
			if (!bdashprune) {														// if branch b not pruned, proceed
				// ----- Generate combos and copy branches
				uvec lmask(cl_layer == Qdash.l);									// layer mask (nodes in current layer)
				vector<uvec> combos;
				uvec bdashvec(1, fill::value(Qdash.b));								// bdash vector, initialise with current branch
				gencombos(combos, Uhat.col(Qdash.b-1), U.col(Qdash.b-1), numnodes, Qdash.k, lmask,roots);

				if (combos.size() > 1) {											// more than 1 branch - create new branches (bmax+1):(bmax+combos.size()-1)
					bdashvec = join_cols(bdashvec, linspace<uvec>(bmax + 1, bmax + combos.size() - 1, combos.size() - 1));
					bmax += combos.size() - 1;
					for (uword i = 1; i < combos.size(); ++i) {// copy branches
						uword bdash = bdashvec(i) - 1;								// for indexing starting at 0
						Uhat.col(bdash) = Uhat.col(Qdash.b-1);
						U.col(bdash) = U.col(Qdash.b - 1);
						G.col(bdash) = G.col(Qdash.b - 1);
						J.col(bdash) = J.col(Qdash.b - 1);
						br_pl(bdash) = -br_pl(Qdash.b - 1);								// copy branch but do not activate it
						//cout << "brpl " << br_pl(Qdash.b-1) << endl;
						Qdash.batchcp_bbrpl(Q, Qdash.b, Qdash.l, bdash+1);			// bdash+1 to get branch id
						Q.batchcp_bbrpl(Q, Qdash.b, Qdash.l, bdash + 1);			// bdash+1 to get branch id
					}
				} // end if copy branches


				uvec indghat(linspace<uvec>(0, combos.size()-1, combos.size()));
				//uword first = 1;
				//cout << "Iter " << iter << " -------------------------------------- " << endl;
				//indghat.print("indghat: ");
				//ghatvec.print("ghat: ");
				//bdashvec.print("bdashvec: ");
				// ---------------- end 202409 ---------------------------------------------------
				//for (auto xind = combos.begin(); xind != combos.end(); ++xind) {	// foreach combo xind
				for (auto ind = indghat.begin(); ind != indghat.end(); ++ind) {	// foreach combo xind
					//vector<vector<uword>> qadd1(numnodes*2, vector<uword>(2));		// qadd: cluster-layer (worse case num parents --> numnodes, sameclust or newclust so *2)
					umat qadd1(numnodes * 2, 2);									// qadd: cluster-layer (worse case num parents --> numnodes, sameclust or newclust so *2)
					uword qaddind = 0;												// index to qadd1
					uvec xind(combos[*ind]);										// xind for that combo
					//xind.print("xind: ");											// debugging
					// ----- Compute G, ghat, update U and J
					if ((xind).size() > 0) {
						uword bdash = bdashvec(*ind); //   test no ghat sorting
						//uword bdash = bdashvec(ind - indghat.begin()); //   with ghat sorting      
						//cout << "b: " << bdash << " brpl: " << br_pl(bdash-1) << endl;
						//if (*ind == 0) {
						//	U.col(bdash - 1).print("U: ");
							//if (Qdash.l > 0)
							//	br_lk.col(Qdash.l - 1 + (bdash - 1) * (lmax + 1)).print("br_lk: ");
						//}
						//if (ind-indghat.begin() == 0)
						//	U.col(bdash - 1).print("U: ");
						uvec intlk(get_intlk(xind, bdash, Qdash.k));				// get intlk mapping
						vec Jdash(numnodes);										// Jdash for propchild and bwdselimcost to update J
						// cost from propchild and bwdselimcost to update G(X_l)
						double cost = propchild(xind, dag, uvec(U.col(bdash - 1)), numnodes, cl_layer, numstates, J.col(bdash - 1), Jdash);
						cost += bwdselimcost((xind)(sort_index(cl_order(xind))), dag, intlk, numstates, Jdash);

						if (Qdash.l > 0 && G(Qdash.l, bdash - 1) == 0)				// if first cluster in layer; add previous layer G
							G(Qdash.l, bdash - 1) = G(Qdash.l - 1, bdash - 1);
						G(Qdash.l, bdash - 1) += cost;								// add clust-layer cost

						//if (G(Qdash.l, bdash - 1) != Gvec(*ind))
						//	cout << "iter " << iter << " not equal, G is " << G(Qdash.l, bdash - 1) << " Gvec is " << Gvec(ind - indghat.begin()) << endl;
						//else
						//	cout << "iter " << iter << " equal, huzzah" << endl;

						setspmscal(U, xind, bdash - 1, Qdash.k);					// set U(xind,bdash)=Qdash.k
						uword indsj = Qdash.l * numnodes * numnodes + (Qdash.k - 1) * numnodes;
						setspmvec(J, linspace<uvec>(indsj, indsj + numnodes - 1, numnodes), bdash - 1, Jdash);// set J

						uvec Uzeros(U.col(bdash - 1));								// heuristic on U(,branch)==0; convert to full vector
						double ghat = heuristic(find(Uzeros == (uword)0), dag, Uzeros, numnodes,
							cl_layer, cl_order, numstates, J.col(bdash - 1)) + G(Qdash.l, bdash - 1);// get heuristic cost
						// Update debugdat
						debugdat(span(0, 6), debugi) = { (double)iter, (double)bdash, (double)Qdash.k, (double)Qdash.l, G(Qdash.l, bdash - 1), ghat, Gmin };
						debugdat(span(8, 7 + numnodes), debugi) = conv_to<vec>::from(Uzeros);
						debugi++;

						uvec lmaski = find(lmask);									// lmask = cl_layer == Qdash.l
						// ----- Layer updates
						if (findspcol_allnz(U.col(bdash - 1), lmaski)) {			// all non-zero U(cl_layer==Qdash.l, bdash-1)
							//cout << "layer " << Qdash.l << endl;
							//Uzeros(lmaski).print("Uzeros: ");
							br_pl(bdash - 1) = (sword)Qdash.l + 1;					// now popped up to layer l
							if (Qdash.l == lmax) {									// termination handling
								if (G(Qdash.l, bdash - 1) < Gmin) {
									Gminind = iter;									// update index of last change in Gmin
									Gmin = G(Qdash.l, bdash - 1);					// update Gmin
									Q.batchrm_ggmin(Gmin);							// flush any dominated Q entries
									Qdash.batchrm_ggmin(Gmin);						// flush any dominated Q entries
								}
							}
							uvec udash(spcolsubset(U.col(bdash - 1), lmaski));		// U of current branch-layer
							for (auto it = br_pl.begin(); it != br_pl.end(); ++it) {// update br_lk
								if (*it >= (sword)(Qdash.l + 1)) {					// if completed up to this layer or later
									uvec u(spcolsubset(U.col(it.row()), lmaski));	// u of branch on br_pl
									if (all(udash == u)) {							// update br_lk (link bdash and b=it.row()+1)
										br_lk(Qdash.l + (bdash - 1) * (lmax + 1), Qdash.l + it.row() * (lmax + 1)) = (uword)1;
										br_lk(Qdash.l + it.row() * (lmax + 1), Qdash.l + (bdash - 1) * (lmax + 1)) = (uword)1;
									}
								} // end if completed up to this layer
							} // end foreach non-zero br_pl
							Q.batchmv_bbrpl(Qdash, bdash, (uword)br_pl(bdash - 1));	// move Q to Qdash any entries activated now due to layer update
							Q.batchrm_bbrpl(bdash, Qdash.l);						// remove Q entries <=l on bdash
							Qdash.batchrm_bbrpl(bdash, Qdash.l);
						} // end layer update

						// ----- Propose parind clusts
						if (G(Qdash.l, bdash - 1) < Gmin) {								// if G < Gmin, proceed
							uvec parinds = unique(mod(find((dag.rows(xind) >= 0).t()), numnodes));// get parents of xind nodes
							for (auto itp = parinds.begin(); itp != parinds.end(); ++itp) {	// for each parind
								uvec scinds(find((sameclust.col(*itp) > 0)));
								uvec parclusts(spcolsubset(U.col(bdash - 1), scinds));// generate sameclust for parclust
								parclusts = join_cols(parclusts, uvec(1, fill::value(cl_order(*itp) + 1)));
								parclusts = GetNondupCol(parclusts);// add to that new clust
								//parinds.print("parinds: ");
								//scinds.print("scinds: ");
								//parclusts.print("parclusts: ");
								uvec indsu((parclusts - 1) * numnodes + *itp);			// update Uhat(indsu,bdash-1)=1
								//indsu.print("indsu: ");
								setspmscal(Uhat, indsu, bdash - 1, (uword)1);
								for (auto itc = parclusts.begin(); itc != parclusts.end(); ++itc) {
									qadd1(qaddind, 0) = *itc;						// col1=clust
									qadd1(qaddind, 1) = cl_layer(*itp);				// col2=layer
									//qadd1[qaddind][0] = *itc;				// col1=clust
									//qadd1[qaddind][1] = cl_layer(*itp);		// col2=layer
									qaddind++;
								}
							} // end foreach parind
						} // end if G < Gmin for parind proposals
						// Add parclust proposals to Q or Qdash - as xind unique to each branch, do separately
						//qadd1.resize(qaddind);									// get rid of zero rows
						//qadd1 = GetUniqueRows(qadd1);								// unique clust-layer proposals, add if not already on Q
						if (qaddind > 0) {
							qadd1 = qadd1.rows(0, qaddind - 1);						// get rid of zero rows
							//qadd1.print("qadd1before: ");
							qadd1 = GetNondupRows(qadd1);							// unique clust-layer proposals, add if not already on Q
							//qadd1.print("qadd1: ");
							Qdash.addtoq(bdash, qadd1, ghat, G(Qdash.l, bdash - 1), (double)br_pl(bdash - 1));	// current branch active for parclust proposal
							Q.addtoq(bdash, qadd1, ghat, G(Qdash.l, bdash - 1), lmax + 1);// current branch inactive for parclust proposal
						}
						//for (uword i = 0; i < qadd1.size(); ++i)					// foreach qadd1 entry
						//	cout << qadd1[i][0] << "," << qadd1[i][1] << endl;
					} // end if xind size > 0
				} // end foreach xind combo
				//} // end for not empty combos and no single element combos[0]={}
			} // end if branch was pruned
		} else {																	// Qdash and Q both empty
			break;
		} // end if
		iter++;

	} // end while
	
	cout << endl << endl << "____________________________________________________" << endl;

	debugdat = debugdat.t();
	debugdat = debugdat.rows(0, debugi - 1);

}




