Rs = QQ[s, u];
Rt = QQ[t, v];
R = Rs ** Rt;


--TSPLINE DEFINITION--
faceDeficits = {1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
edgeDeficits = {1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
vertexDeficits = {1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
smoothness = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
maxDegree = 20;


--COMPLEX INTERSECTION--
INT = (myComplex, j) -> intersect(myComplex, directSum apply(numgens ambient myComplex, i-> module ideal ((u*v)^j)));
INTMOD = (myComplex, j) -> ((INT(myComplex,j))/(INT(myComplex,j+1)));
EULER = (myComplex0, myComplex1, myComplex2, j) -> hilbertFunction({j, j}, myComplex2) - hilbertFunction({j, j}, myComplex1) + hilbertFunction({j, j}, myComplex0);


--BOUNDARY MAPS--
D2 = matrix(R,
			{{1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0},
			{0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1}});
D1 = matrix(R,
			{{1, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 1, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 1, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, -1},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1}});


--CONSTANT COMPLEX--
C2 = directSum apply(#faceDeficits, i -> module ideal((u*v)^(faceDeficits_i)));
C1 = directSum apply(#edgeDeficits, i -> module ideal((u*v)^(edgeDeficits_i)));
C0 = directSum apply(#vertexDeficits, i -> module ideal((u*v)^(vertexDeficits_i)));
C = chainComplex(inducedMap(C0, C1, D1), inducedMap(C1, C2, D2));
C2'_0 = INT(C2, 1);C1'_0 = INT(C1, 1);C0'_0 = INT(C0, 1);C'_0 = chainComplex(inducedMap(C0'_0, C1'_0, D1), inducedMap(C1'_0, C2'_0, D2));C2'_1 = INTMOD(C2, 0);C1'_1 = INTMOD(C1, 0);C0'_1 = INTMOD(C0, 0);C'_1 = chainComplex(inducedMap(C0'_1, C1'_1, D1), inducedMap(C1'_1, C2'_1, D2));

--IDEALS COMPLEX--
linearForms = {
		(s-u*1/4)^(3),
		(s-u*1/2)^(3),
		(s-u*3/4)^(3),
		(s-u*1/4)^(3),
		(t-v*1/4)^(3),
		(s-u*1/2)^(3),
		(t-v*1/4)^(3),
		(s-u*3/4)^(3),
		(t-v*1/4)^(3),
		(t-v*1/4)^(3),
		(s-u*1/4)^(3),
		(t-v*3/8)^(3),
		(s-u*1/2)^(3),
		(t-v*3/8)^(3),
		(s-u*3/4)^(3),
		(t-v*3/8)^(3),
		(t-v*3/8)^(3),
		(s-u*1/4)^(3),
		(t-v*1/2)^(3),
		(s-u*1/2)^(3),
		(t-v*1/2)^(3),
		(s-u*3/4)^(3),
		(t-v*1/2)^(3),
		(t-v*1/2)^(3),
		(s-u*1/4)^(3),
		(t-v*5/8)^(3),
		(s-u*1/2)^(3),
		(t-v*5/8)^(3),
		(s-u*3/4)^(3),
		(t-v*5/8)^(3),
		(t-v*5/8)^(3),
		(s-u*1/4)^(3),
		(t-v*3/4)^(3),
		(s-u*1/2)^(3),
		(t-v*3/4)^(3),
		(s-u*3/4)^(3),
		(t-v*3/4)^(3),
		(t-v*3/4)^(3),
		(t-v*1/4)^(3),
		(t-v*3/8)^(3),
		(s-u*3/8)^(3),
		(t-v*1/2)^(3),
		(s-u*3/8)^(3),
		(t-v*5/8)^(3),
		(s-u*3/8)^(3),
		(t-v*3/4)^(3),
		(s-u*3/8)^(3),
		(t-v*1/4)^(3),
		(t-v*3/8)^(3),
		(s-u*5/8)^(3),
		(t-v*1/2)^(3),
		(s-u*5/8)^(3),
		(t-v*5/8)^(3),
		(s-u*5/8)^(3),
		(t-v*3/4)^(3),
		(s-u*5/8)^(3),
		(t-v*1/4)^(3),
		(t-v*3/8)^(3),
		(s-u*7/8)^(3),
		(t-v*1/2)^(3),
		(s-u*7/8)^(3),
		(t-v*5/8)^(3),
		(s-u*7/8)^(3),
		(t-v*3/4)^(3),
		(s-u*7/8)^(3)};
I2 = image matrix{{0_R}};
I1 = directSum apply(#edgeDeficits, i -> module ideal((linearForms_i)*((u*v)^(edgeDeficits_i))));
I0 = directSum {
		module ideal((linearForms_0)*((u*v)^(edgeDeficits_0)), (linearForms_3)*((u*v)^(edgeDeficits_3)), (linearForms_4)*((u*v)^(edgeDeficits_4)), (linearForms_6)*((u*v)^(edgeDeficits_6))),
		module ideal((linearForms_1)*((u*v)^(edgeDeficits_1)), (linearForms_5)*((u*v)^(edgeDeficits_5)), (linearForms_8)*((u*v)^(edgeDeficits_8)), (linearForms_38)*((u*v)^(edgeDeficits_38))),
		module ideal((linearForms_2)*((u*v)^(edgeDeficits_2)), (linearForms_7)*((u*v)^(edgeDeficits_7)), (linearForms_9)*((u*v)^(edgeDeficits_9)), (linearForms_47)*((u*v)^(edgeDeficits_47))),
		module ideal((linearForms_3)*((u*v)^(edgeDeficits_3)), (linearForms_10)*((u*v)^(edgeDeficits_10)), (linearForms_11)*((u*v)^(edgeDeficits_11)), (linearForms_13)*((u*v)^(edgeDeficits_13))),
		module ideal((linearForms_5)*((u*v)^(edgeDeficits_5)), (linearForms_12)*((u*v)^(edgeDeficits_12)), (linearForms_15)*((u*v)^(edgeDeficits_15)), (linearForms_39)*((u*v)^(edgeDeficits_39))),
		module ideal((linearForms_7)*((u*v)^(edgeDeficits_7)), (linearForms_14)*((u*v)^(edgeDeficits_14)), (linearForms_16)*((u*v)^(edgeDeficits_16)), (linearForms_48)*((u*v)^(edgeDeficits_48))),
		module ideal((linearForms_10)*((u*v)^(edgeDeficits_10)), (linearForms_17)*((u*v)^(edgeDeficits_17)), (linearForms_18)*((u*v)^(edgeDeficits_18)), (linearForms_20)*((u*v)^(edgeDeficits_20))),
		module ideal((linearForms_12)*((u*v)^(edgeDeficits_12)), (linearForms_19)*((u*v)^(edgeDeficits_19)), (linearForms_22)*((u*v)^(edgeDeficits_22)), (linearForms_41)*((u*v)^(edgeDeficits_41))),
		module ideal((linearForms_14)*((u*v)^(edgeDeficits_14)), (linearForms_21)*((u*v)^(edgeDeficits_21)), (linearForms_23)*((u*v)^(edgeDeficits_23)), (linearForms_50)*((u*v)^(edgeDeficits_50))),
		module ideal((linearForms_17)*((u*v)^(edgeDeficits_17)), (linearForms_24)*((u*v)^(edgeDeficits_24)), (linearForms_25)*((u*v)^(edgeDeficits_25)), (linearForms_27)*((u*v)^(edgeDeficits_27))),
		module ideal((linearForms_19)*((u*v)^(edgeDeficits_19)), (linearForms_26)*((u*v)^(edgeDeficits_26)), (linearForms_29)*((u*v)^(edgeDeficits_29)), (linearForms_43)*((u*v)^(edgeDeficits_43))),
		module ideal((linearForms_21)*((u*v)^(edgeDeficits_21)), (linearForms_28)*((u*v)^(edgeDeficits_28)), (linearForms_30)*((u*v)^(edgeDeficits_30)), (linearForms_52)*((u*v)^(edgeDeficits_52))),
		module ideal((linearForms_24)*((u*v)^(edgeDeficits_24)), (linearForms_31)*((u*v)^(edgeDeficits_31)), (linearForms_32)*((u*v)^(edgeDeficits_32)), (linearForms_34)*((u*v)^(edgeDeficits_34))),
		module ideal((linearForms_26)*((u*v)^(edgeDeficits_26)), (linearForms_33)*((u*v)^(edgeDeficits_33)), (linearForms_36)*((u*v)^(edgeDeficits_36)), (linearForms_45)*((u*v)^(edgeDeficits_45))),
		module ideal((linearForms_28)*((u*v)^(edgeDeficits_28)), (linearForms_35)*((u*v)^(edgeDeficits_35)), (linearForms_37)*((u*v)^(edgeDeficits_37)), (linearForms_54)*((u*v)^(edgeDeficits_54))),
		module ideal((linearForms_6)*((u*v)^(edgeDeficits_6)), (linearForms_38)*((u*v)^(edgeDeficits_38)), (linearForms_40)*((u*v)^(edgeDeficits_40))),
		module ideal((linearForms_13)*((u*v)^(edgeDeficits_13)), (linearForms_39)*((u*v)^(edgeDeficits_39)), (linearForms_40)*((u*v)^(edgeDeficits_40)), (linearForms_42)*((u*v)^(edgeDeficits_42))),
		module ideal((linearForms_20)*((u*v)^(edgeDeficits_20)), (linearForms_41)*((u*v)^(edgeDeficits_41)), (linearForms_42)*((u*v)^(edgeDeficits_42)), (linearForms_44)*((u*v)^(edgeDeficits_44))),
		module ideal((linearForms_27)*((u*v)^(edgeDeficits_27)), (linearForms_43)*((u*v)^(edgeDeficits_43)), (linearForms_44)*((u*v)^(edgeDeficits_44)), (linearForms_46)*((u*v)^(edgeDeficits_46))),
		module ideal((linearForms_34)*((u*v)^(edgeDeficits_34)), (linearForms_45)*((u*v)^(edgeDeficits_45)), (linearForms_46)*((u*v)^(edgeDeficits_46))),
		module ideal((linearForms_8)*((u*v)^(edgeDeficits_8)), (linearForms_47)*((u*v)^(edgeDeficits_47)), (linearForms_49)*((u*v)^(edgeDeficits_49))),
		module ideal((linearForms_15)*((u*v)^(edgeDeficits_15)), (linearForms_48)*((u*v)^(edgeDeficits_48)), (linearForms_49)*((u*v)^(edgeDeficits_49)), (linearForms_51)*((u*v)^(edgeDeficits_51))),
		module ideal((linearForms_22)*((u*v)^(edgeDeficits_22)), (linearForms_50)*((u*v)^(edgeDeficits_50)), (linearForms_51)*((u*v)^(edgeDeficits_51)), (linearForms_53)*((u*v)^(edgeDeficits_53))),
		module ideal((linearForms_29)*((u*v)^(edgeDeficits_29)), (linearForms_52)*((u*v)^(edgeDeficits_52)), (linearForms_53)*((u*v)^(edgeDeficits_53)), (linearForms_55)*((u*v)^(edgeDeficits_55))),
		module ideal((linearForms_36)*((u*v)^(edgeDeficits_36)), (linearForms_54)*((u*v)^(edgeDeficits_54)), (linearForms_55)*((u*v)^(edgeDeficits_55))),
		module ideal((linearForms_9)*((u*v)^(edgeDeficits_9)), (linearForms_56)*((u*v)^(edgeDeficits_56)), (linearForms_58)*((u*v)^(edgeDeficits_58))),
		module ideal((linearForms_16)*((u*v)^(edgeDeficits_16)), (linearForms_57)*((u*v)^(edgeDeficits_57)), (linearForms_58)*((u*v)^(edgeDeficits_58)), (linearForms_60)*((u*v)^(edgeDeficits_60))),
		module ideal((linearForms_23)*((u*v)^(edgeDeficits_23)), (linearForms_59)*((u*v)^(edgeDeficits_59)), (linearForms_60)*((u*v)^(edgeDeficits_60)), (linearForms_62)*((u*v)^(edgeDeficits_62))),
		module ideal((linearForms_30)*((u*v)^(edgeDeficits_30)), (linearForms_61)*((u*v)^(edgeDeficits_61)), (linearForms_62)*((u*v)^(edgeDeficits_62)), (linearForms_64)*((u*v)^(edgeDeficits_64))),
		module ideal((linearForms_37)*((u*v)^(edgeDeficits_37)), (linearForms_63)*((u*v)^(edgeDeficits_63)), (linearForms_64)*((u*v)^(edgeDeficits_64)))};
I = chainComplex(inducedMap(I0, I1, D1));
I1'_0 = INT(I1, 1);I0'_0 = INT(I0, 1);I'_0 = chainComplex(inducedMap(I0'_0, I1'_0, D1));I1'_1 = INTMOD(I1, 0);I0'_1 = INTMOD(I0, 0);I'_1 = chainComplex(inducedMap(I0'_1, I1'_1, D1));

--QUOTIENT COMPLEX--
Q2 = C2;
Q1 = C1 / I1;
Q0 = C0 / I0;
Q = chainComplex(inducedMap(Q0, Q1, D1), inducedMap(Q1, Q2, D2));
Q2'_0 = INT(C2, 1);Q1'_0 = INT(C1, 1) / INT(I1, 1);Q0'_0 = INT(C0, 1) / INT(I0, 1);Q'_0 = chainComplex(inducedMap(Q0'_0, Q1'_0, D1), inducedMap(Q1'_0, Q2'_0, D2));Q2'_1 = INTMOD(C2, 0);Q1'_1 = INTMOD(C1, 0) / INTMOD(I1, 0);Q0'_1 = INTMOD(C0, 0) / INTMOD(I0, 0);Q'_1 = chainComplex(inducedMap(Q0'_1, Q1'_1, D1), inducedMap(Q1'_1, Q2'_1, D2));

--DIMENSIONS OF HOMOLOGIES--
--Constant
H = homology(C);
c2 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_2) );
c1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_1) );
c0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_0) );
H = homology(C'_0);
c2'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_2) );
c1'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_1) );
c0'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_0) );
H = homology(C'_1);
c2'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_2) );
c1'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_1) );
c0'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_0) );


--Ideal
H = homology(I);
i1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_1) );
i0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_0) );
H = homology(I'_0);
i1'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_1) );
i0'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_0) );
H = homology(I'_1);
i1'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_1) );
i0'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_0) );


--Quotient
H = homology(Q);
SpM2 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_2) );
SpM2Int'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, INT(H_2, 1)) );
SpM2Int'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, INTMOD(H_2, 0)) );
q1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_1) );
q0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_0) );
qEuler = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, Q2) - hilbertFunction({i, i}, Q1) + hilbertFunction({i, i}, Q0));
SpEuler = apply(toList(1..maxDegree), i-> EULER(Q0, Q1, Q2, i) );
SpEuler = SpEuler + i0 - c0;
H = homology(Q'_0);
SpM2Deg'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_2) );
q1'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_1) );
q0'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_0) );
q00'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, Q0'_0) );
q11'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, Q1'_0) );
q22'_0 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, Q2'_0) );
H = homology(Q'_1);
SpM2Deg'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_2) );
q1'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_1) );
q0'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, H_0) );
q00'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, Q0'_1) );
q11'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, Q1'_1) );
q22'_1 = apply(toList(1..maxDegree), i-> hilbertFunction({i, i}, Q2'_1) );


--Output
dl = prepend("Degree ", toList(1..maxDegree));
qEuler = prepend("Euler(Q) [M2] ", qEuler);
q1 = prepend("H_1(Q) [M2] ", q1);
q0 = prepend("H_0(Q) [M2] ", q0);
q1'_0 = prepend("H_1(Q)'_0 [M2] ", q1'_0);
q0'_0 = prepend("H_0(Q)'_0 [M2] ", q0'_0);
q00'_0 = prepend("Q0'_0 [M2] ", q00'_0);
q11'_0 = prepend("Q1'_0 [M2] ", q11'_0);
q22'_0 = prepend("Q2'_0 [M2] ", q22'_0);
q1'_1 = prepend("H_1(Q)'_1 [M2] ", q1'_1);
q0'_1 = prepend("H_0(Q)'_1 [M2] ", q0'_1);
q00'_1 = prepend("Q0'_1 [M2] ", q00'_1);
q11'_1 = prepend("Q1'_1 [M2] ", q11'_1);
q22'_1 = prepend("Q2'_1 [M2] ", q22'_1);
c2 = prepend("H_2(C) [M2] ", c2);
c1 = prepend("H_1(C) [M2] ", c1);
c0 = prepend("H_0(C) [M2] ", c0);
c2'_0 = prepend("H_2(C)'_0 [M2] ", c2'_0);
c1'_0 = prepend("H_1(C)'_0 [M2] ", c1'_0);
c0'_0 = prepend("H_0(C)'_0 [M2] ", c0'_0);
c2'_1 = prepend("H_2(C)'_1 [M2] ", c2'_1);
c1'_1 = prepend("H_1(C)'_1 [M2] ", c1'_1);
c0'_1 = prepend("H_0(C)'_1 [M2] ", c0'_1);
i1 = prepend("H_1(I) [M2] ", i1);
i0 = prepend("H_0(I) [M2] ", i0);
i1'_0 = prepend("H_1(I)'_0 [M2] ", i1'_0);
i0'_0 = prepend("H_0(I)'_0 [M2] ", i0'_0);
i1'_1 = prepend("H_1(I)'_1 [M2] ", i1'_1);
i0'_1 = prepend("H_0(I)'_1 [M2] ", i0'_1);
SpM2 = prepend("Sp [M2] ", SpM2);

SpEuler = prepend("Sp [Euler] ", SpEuler);

SpDeg = prepend("SpDeg [M2] ", SpM2Deg'_0 + SpM2Deg'_1);

SpM2Deg'_0 = prepend("S_0 [M2] ", SpM2Deg'_0);

SpM2Deg'_1 = prepend("S_1 [M2] ", SpM2Deg'_1);

netList {dl, qEuler, q1, q0}
netList {dl, c2, c1, c0}
netList {dl, c2'_0, c1'_0, c0'_0}
netList {dl, c2'_1, c1'_1, c0'_1}
netList {dl, i1, i0}
netList {dl, i1'_0, i0'_0}
netList {dl, i1'_1, i0'_1}
netList {dl, c0'_0, i0'_0}
netList {dl, c0'_1, i0'_1}
netList {dl, i0'_0, c0'_0, q1'_0, q0'_0}
netList {dl, i0'_1, c0'_1, q1'_1, q0'_1}
netList {dl, SpM2Deg'_0, SpM2Deg'_1}
netList {dl, q22'_0, q22'_1}
netList {dl, q11'_0, q11'_1}
netList {dl, q00'_0, q00'_1}
netList {dl, SpM2, SpEuler, SpDeg}
