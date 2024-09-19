nthreads=17

for((i=1; i<=$nthreads; i++))
do	
	output_file=out_thread_"$i".txt
	magma  thread_id:=$i exp_various_l_finding_vertices_in_Fp_S_1.m > $output_file &
done


