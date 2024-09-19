p:= 35184372088777
;
terms_irr_poly:= [5, 0, 1]
;
seq_j_inv:= [[32958660381388, 29001661615062], [34894147070814, 8997070700500], [29900321453667, 15909095809544], [18973985527669, 14731231721230], [24929578619126, 22325974622728], [30635091681, 24762878180618], [30041548431956, 21331880121464], [19193371496766, 30230577030923], [29398540395252, 34008196813121], [24365201634513, 2528621119893], [21277385847765, 29148928438654], [811711818310, 12078500168894], [24438569884588, 14906951612472], [21075941728705, 7223610757145], [6443574753755, 32648085111492], [527677967013, 13679662440141], [22773658472974, 32418891523550], [12866521962465, 5029563964572], [13649594266392, 10106383509254], [25236632126188, 7406938422955]]
;
Number_of_curves:=#seq_j_inv;
Fp:=GF(p);
Fp2:=ExtensionField<Fp,z|&+[terms_irr_poly[i]*z^(i-1): i in [1..#terms_irr_poly]]>;

printf "\n p:=%o;",p;
printf "\n Defining polynomial=%o",DefiningPolynomial(Fp2);


R2<x,y>:=PolynomialRing(Fp2,2);


l:=PrimesInInterval(2,59)[StringToInteger(thread_id)];

mod_poly:=R2!ClassicalModularPolynomial(l);//The modular polynomial construction
never_walk_j_inv:=99*p^2+99;
printf "\n time to construct %o cycles FOR DIST=2 method and l=%o is",Number_of_curves,l;



time for i in [1..#seq_j_inv] do
	j0:=seq_j_inv[i][1]+seq_j_inv[i][2]*Fp2.1;	
	j:=j0;
	signal:=0;
	j_previous:=never_walk_j_inv;
	while(signal eq 0) do
		f:=Evaluate(mod_poly,y,j);
		set_of_roots:=Roots(UnivariatePolynomial(f));
		all_roots:=[set_of_roots[i][1] : i in [1..#set_of_roots]];
		Exclude(~all_roots,j_previous);
		for t in [1..#all_roots] do
			j1:=all_roots[t];		
			j1_power_p:=j1^p;
			if(j1 eq j1_power_p) then
				break;
			else
				event1:=Evaluate(mod_poly,[j1,j1_power_p]) eq 0;	
				if(event1) then
					signal:=1;
					break;
				else
					f1:=Evaluate(mod_poly,y,j1);
					f2:=Evaluate(mod_poly,y,j1_power_p);
				event2:=GreatestCommonDivisor(UnivariatePolynomial(f1),UnivariatePolynomial(f2)) ne 1;
					if(event2) then
						signal:=1;
						break;
					end if;	
				end if;
			end if;
		end for;
		j_previous:=j;		
		j:=Random(all_roots);
	end while;	
end for;	
