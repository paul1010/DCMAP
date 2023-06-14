// dcmap_cost.h stores cost and related functions
#define MARGFRAC 0.6


// ----- propchild
// Propagate J from child clust-layers to nodes in xind
// Inputs: uvec xind = indices of nodes[1:n] to propagate to, dag = cptparchild,
// U = uvec cluster mapping for this branch U[,b], n=numnodes, cl.layer, numstates foreach node
// J = vec [k,l,i] for this branch J[,b]; vec of length n*lmax*n (joint distribution)
// Outputs: Jdash = updated joint distribution + returns cost
double propchild(const uvec& xind, const mat& dag, const uvec& U, uword n, const uvec& cl_layer, const uvec& numstates, 
	SpSubview_col<double>& J, vec& Jdash) {
	Jdash.zeros();												// joint distribution and mask version of x inds
	vec Jcpt(n);												// joint of CPTs of nodes in xind (mask)
	vec Zero(n, fill::value(0));								// zero vec for logical comparisons with Jcpt, childJ
	double cost = 0;											// computational cost
	uvec xc;													// child nodes of xind (as index)

	Jcpt(unique(join_cols(unique(mod(find((dag.rows(xind) >= 0).t()), n)), xind))).fill(1);	// Joint CPT nodes in xind
	xc = GetNondupVec(mod(find((dag.cols(xind) >= 0)), n));
	xc = std_setdiff(xc, xind); //						       
	//xc.print("xc: ");										//debugging
	//Jcpt(2) = 0;												// debugging
	//Jcpt(1) = 0;												// debugging
	////Jcpt(0) = 0;											// debugging
	//Jcpt.print("Jcpt: ");										//debugging
	if (xc.size() > 0) {										// if there are child clust-layers to propagate
		//vector<vector<uword>> xc_lk(xc.size(), vector<uword>(2)); // get unique clust-layers from xc
		//for (uword i = 0; i < xc.size(); i++) {
		//	xc_lk[i][0] = cl_layer(xc(i));						// col0=layer
		//	xc_lk[i][1] = U(xc(i));								// col1=clust
		//}
		//xc_lk = GetUniqueRows(xc_lk);							// unique clust-layers
		
		umat xc_lk(xc.size(), 2);								// get unique clust-layers from xc (non-duplicated to mimic R unique())
		xc_lk.col(0) = cl_layer(xc);							// col0=layer
		xc_lk.col(1) = U(xc);									// col1=clust
		xc_lk = GetNondupRows(xc_lk);
		//xc_lk.print("xc_lk: ");

		//J(2) = 1;												//debugging
		//J(42) = 1;											//debugging
		//J(43) = 1;											//debugging
		//xc_lk[0][1] = 1;										//debugging
		for (uword i = 0; i < xc_lk.n_rows; i++) {	// foreach uniq child clust-layer xc_lk.size() for vector<> implementation
			// Get child clust-layer J
			//uword indsj = xc_lk[i][0] * n * n + (xc_lk[i][1] - 1) * n;
			uword indsj = xc_lk(i,0) * n * n + (xc_lk(i,1) - 1) * n;
			//cout << xc_lk[i][0] << "," << xc_lk[i][1] << endl;//debugging
			//cout << "indsj = " << indsj << endl;				//debugging
			vec childJ(J.rows(indsj, indsj + n - 1));
			//childJ.print("childJ: ");							// debugging

			// Marg out all nodes in childJ not in CPTs(xind)
			uvec margnodesind = find(childJ != Zero && Jcpt == Zero);
			uvec multnodesind = find(childJ != Zero && Jcpt != Zero);
			// by de Morgan's rule !(A&B)=!A|!B, so childJ & (!childJ || Jcpt) = childJ & !childJ || childJ & Jcpt = childJ & Jcpt
			//margnodesind.print("margnodesind: ");				// debugging
			//multnodesind.print("multnodesind: ");
			double jbase = prod(numstates(multnodesind));
			// cumulative product of numstates(margnodesmask) w. last node = numstates-1
			umat NS(margnodesind.size(), margnodesind.size());
			NS.each_col() = numstates(margnodesind);
			NS.diag() -= 1;
			if (margnodesind.size() > 1)
				NS.elem(trimatl_ind(size(NS), -1)).ones();		// change lower triangular without diag into 1
			cost += (double)sum(prod(NS, 0)) * jbase * MARGFRAC;// compute cost

			//Jdash.print("Jdash: ");							// debugging
			//cout << "jbase: " << jbase << endl;				// debugging
			//NS.print("NS: ");									// debugging
			//cout << "cost: " << cost << endl;					// debugging

			// Multiply
			Jdash(multnodesind).ones();
			cost += prod(numstates(find(Jdash)));
			//cout << "cost mult: " << cost << endl;
		}
	}
	return cost;
}


// ----- bwdselimcost
// Multiply and marginalise nodes, finding cost and effect on J
// Inputs: uvec xind = indices of nodes[1:n] to do bwds elimination in order to mult-marg (eliminate), 
// dag = cptparchild, uvec intlk(xind) 1=int,2=lk node
// uvec numstates foreach node, Jdash initial J distribution (size numnodes) updated at end of mult-marg
// Outputs: updated joint distribution Jdash, and cost
double bwdselimcost(const uvec& xind, const mat& dag, const uvec& intlk, const uvec& numstates, vec& Jdash) {
	double cost = 0;
	for (uword i = 0; i < xind.size(); i++) {
		// multiply
		Jdash(find((dag.row(xind(i)) >= 0))).ones();	// add node parents
		Jdash(xind(i)) = 1.0;							// add node itself
		cost += prod(numstates(find(Jdash > 0)));
		// marginalise
		if (intlk(i) == 1) {
			uvec ns(numstates);
			ns(xind(i)) -= 1;
			cost += MARGFRAC * prod(ns(find(Jdash > 0)));
			Jdash(xind(i)) = 0;
		}
		//Jdash.print("Jdash: ");
	}
	return cost;
}


// ----- heuristic
// Compute heuristic estimate of cost given nodes and current J
// Inputs: uvec xind = indices of nodes[1:n] to calc heuristic, dag = cptparchild,
// U = uvec cluster mapping for this branch U[,b], n=numnodes, cl.layer, cl.order, numstates foreach node
// J = vec [k,l,i] for this branch J[,b]; vec of length n*lmax*n (joint distribution)
// Outputs: computational cost estimate as a scalar
double heuristic(const uvec& xind, const mat& dag, const uvec& U, uword n, const uvec& cl_layer, const uvec& cl_order, 
	const uvec& numstates, SpSubview_col<double>& J) {
	//uvec layers = sort(unique(cl_layer(xind)));
	
	double cost = 0;
	vec Jdash(n);
	cost = propchild(xind, dag, U, n, cl_layer, numstates, J, Jdash);
	//cout << "propcost = " << cost << endl;
	//J.print("J: ");
	uvec intlk(xind.size(), fill::ones);
	cost += bwdselimcost(xind(sort_index(cl_order(xind))), dag, intlk, numstates, Jdash);
	return cost;
}

// ----- gencombos
// generate combinations with root node heuristic (no proposals with >1 root node in cluster)
// given Uhat=Uhat[,b], U=U[,b], numnodes, cluster k, uvec lmask of nodes in current layer
// roots uvec contains indices to root nodes for root node heuristic (no 2 root nodes in same cluster)
// combos is the vecotr<uvec> to write xind after generating combos
void gencombos(vector<uvec>& combos, const SpSubview_col<uword>& Uhat, const SpSubview_col<uword>& U, 
		uword numnodes, uword k, const uvec& lmask, const uvec& roots) {
	//vector<vector<uword>> combos;
	//vector<uvec> combos;
	umat uhat(Uhat);												// cluster proposals k-i for branch b as mat
	uhat.reshape(numnodes, numnodes);								// columns are clusters, rows = nodes, 1=member of cluster
	uvec u(U);														// current cluster assignments

	//uhat.print("uhat: ");
	uvec z1(sum(uhat, 1)<=1);										// nodes with no other proposals
	//z1.print("z1: ");
	//lmask.print("lmask: ");;
	
	uvec z(uhat.col(k - 1) && lmask);								// nodes in this cluster-layer
	//z.print("z: ");
	// nodes for cluster-layer that have options/combos for assignment, currently (via U) in k or unassigned
	uvec z2(z && (z1 == 0) && (u == k || u == 0));					// by de Morgan's laws: not (A and B) = (not A) or (not B)
	//z2.print("z2: ");

	uvec x(find(z2 != 0));
	uvec y(find(z1 && z));
	//uvec xx = { 0,1,4 };
	//x = join_cols(x,xx);
	//x.print("x: ");
	//y.print("y: ");

	if (x.size() > 1) {												// generate combos
		combos.resize(pow(2, x.size()));							// resize combos 2^x.size 
		uword ic = 0;												// to index combos
		for (uword i = 1; i <= x.size(); ++i) {						// for combos of length 1 to length x.size()
			vector<bool> v(x.size());
			std::fill(v.begin(), v.begin()+i, true);				// use std::prev_permutation to generate combos, based on true/false mask
			do {
				uvec xij(i);										// combo is of length i
				uword k = 0;
				for (int j = 0; j < x.size(); ++j) {				// obtain combo from true/false mask
					if (v[j]) {
						//cout << (j + 1) << " ";
						xij(k) = x(j);
						k++;
					}
				}
				//cout << endl;
				//xij.print("xij: ");
				combos[ic] = join_cols(xij, y);
				//combos[ic].print("combosic: ");
				ic++;
			} while (prev_permutation(v.begin(), v.end()));
		}
		//cout << " -------" << endl;
		//for (uword i = 1; i < pow(2, x.size()); ++i) {
		//	uvec mask(x.size());									// mask = vec of inds to include in this combo
		//	uword ii = i;											// copy of i to shift right
		//	for (uword j = 0; j < x.size(); ++j) {					
		//		mask(j) = ii & 1;									// add to mask(j) the rightmost bit
		//		ii >>= 1;											// right shift by one
		//	}
		//	combos[i-1] = join_cols(x(find(mask != 0)), y);			// combos of x
		//	//combos[i - 1].print("combosi: ");
		//}
		combos[combos.size() - 1] = y;								// append y
	} else {														// no combos
		combos.resize(x.size() + 1);
		for (uword i = 0; i < x.size(); ++i) { 
			combos[i] = uvec(1, fill::value(x(i)));					// individual nodes as a "combo"
		}
		combos[combos.size()-1] = y;								// append y
	}

	 // Displaying the 2D vector
	//cout << "combos: " << endl;
	//for (int i = 0; i < combos.size(); i++) {
	//	combos[i].print();
	//}

	// Filter combos which have more than one root node (via intersect function)
	vector<uvec> combos2;
	std::copy_if(combos.begin(), combos.end(), std::back_inserter(combos2), [&](uvec iv) {return uvec(intersect(roots, iv)).size() <= 1; });
	combos = combos2;

	//// Filter combos which have more than one root node (via intersect function)
	//// https://stackoverflow.com/questions/30611584/how-to-efficiently-delete-elements-from-vector-c
	//auto it = combos.begin();
	//while (it != combos.end()) {
	//	if (uvec(intersect(roots, *it)).size() > 1) {
	//		auto lastit = combos.end() - 1;
	//		if (it != lastit) {
	//			*it = move(*lastit);
	//		} else {
	//			combos.pop_back();
	//			break;
	//		}
	//		combos.pop_back();
	//		//cout << "." << combos.size();
	//		//cout << it - combos.begin();
	//	} else {
	//		//cout << "-";
	//		//cout << it - combos.begin();
	//		++it;
	//	}
	//}
	//cout << endl;
}

//// ----- gencombos
//// generate combinations with root node heuristic (no proposals with >1 root node in cluster)
//// given Uhat=Uhat[,b], U=U[,b], numnodes, cluster k, uvec lmask of nodes in current layer
//// roots uvec contains indices to root nodes for root node heuristic (no 2 root nodes in same cluster)
//// combos is the vecotr<uvec> to write xind after generating combos
//void gencombos(vector<uvec>& combos, const SpSubview_col<uword>& Uhat, const SpSubview_col<uword>& U,
//	uword numnodes, uword k, const uvec& lmask, const uvec& roots) {
//	//vector<vector<uword>> combos;
//	//vector<uvec> combos;
//	umat uhat(Uhat);												// cluster proposals k-i for branch b as mat
//	uhat.reshape(numnodes, numnodes);								// columns are clusters, rows = nodes, 1=member of cluster
//	uvec u(U);														// current cluster assignments
//
//	//uhat.print("uhat: ");
//	uvec z1(sum(uhat, 1) <= 1);										// nodes with no other proposals
//	//z1.print("z1: ");
//	//lmask.print("lmask: ");;
//
//	uvec z(uhat.col(k - 1) && lmask);								// nodes in this cluster-layer
//	//z.print("z: ");
//	// nodes for cluster-layer that have options/combos for assignment, currently (via U) in k or unassigned
//	uvec z2(z && (z1 == 0) && (u == k || u == 0));					// by de Morgan's laws: not (A and B) = (not A) or (not B)
//	//z2.print("z2: ");
//
//	uvec x(find(z2 != 0));
//	uvec y(find(z1 && z));
//	//uvec xx = { 0,1,4 };
//	//x = join_cols(x,xx);
//	//x.print("x: ");
//	//y.print("y: ");
//
//	if (x.size() > 1) {												// generate combos
//		combos.resize(pow(2, x.size()));							// resize combos 2^x.size + 1
//		for (uword i = 1; i < pow(2, x.size()); ++i) {
//			uvec mask(x.size());									// mask = vec of inds to include in this combo
//			uword ii = i;											// copy of i to shift right
//			for (uword j = 0; j < x.size(); ++j) {
//				mask(j) = ii & 1;									// add to mask(j) the rightmost bit
//				ii >>= 1;											// right shift by one
//			}
//			combos[i - 1] = join_cols(x(find(mask != 0)), y);			// combos of x
//		}
//		combos[combos.size() - 1] = y;								// append y
//	}
//	else {														// no combos
//		combos.resize(x.size() + 1);
//		for (uword i = 0; i < x.size(); ++i) {
//			combos[i] = uvec(1, fill::value(x(i)));					// individual nodes as a "combo"
//		}
//		combos[combos.size() - 1] = y;								// append y
//	}
//
//	// Displaying the 2D vector
//   //cout << "combos: " << endl;
//   //for (int i = 0; i < combos.size(); i++) {
//   //	for (int j = 0; j < combos[i].size(); j++)
//   //		cout << combos[i][j] << " ";
//   //	cout << endl;
//   //}
//
//   // Filter combos which have more than one root node (via intersect function)
//   // https://stackoverflow.com/questions/30611584/how-to-efficiently-delete-elements-from-vector-c
//	auto it = combos.begin();
//	while (it != combos.end()) {
//		if (uvec(intersect(roots, *it)).size() > 1) {
//			auto lastit = combos.end() - 1;
//			if (it != lastit) {
//				*it = move(*lastit);
//			}
//			else {
//				combos.pop_back();
//				break;
//			}
//			combos.pop_back();
//			//cout << "." << combos.size();
//			//cout << it - combos.begin();
//		}
//		else {
//			//cout << "-";
//			//cout << it - combos.begin();
//			++it;
//		}
//	}
//	//cout << endl;
//}





// ----- 
//// https://www.geeksforgeeks.org/print-all-possible-combinations-of-r-elements-in-a-given-array-of-size-n/
///* arr[] ---> Input Array
//data[] ---> Temporary array to
//store current combination
//start & end ---> Starting and
//Ending indexes in arr[]
//index ---> Current index in data[]
//r ---> Size of a combination to be printed */
//void ncombr(int arr[], int data[],
//	int start, int end,
//	int index, int r)
//{
//	// Current combination is ready
//	// to be printed, print it
//	if (index == r)
//	{
//		for (int j = 0; j < r; j++)
//			cout << data[j] << " ";
//		cout << endl;
//		return;
//	}
//
//	// replace index with all possible
//	// elements. The condition "end-i+1 >= r-index"
//	// makes sure that including one element
//	// at index will make a combination with
//	// remaining elements at remaining positions
//	for (int i = start; i <= end &&
//		end - i + 1 >= r - index; i++)
//	{
//		data[index] = arr[i];
//		ncombr(arr, data, i + 1,
//			end, index + 1, r);
//	}
//}