#ifndef INC_REQ_newGA
#define INC_REQ_newGA
#include "newGA.hh"
#include <math.h>

skeleton newGA{

	// Problem ---------------------------------------------------------------

	Problem::Problem ():_dimension(0),_matrizCostos(NULL)
	{}

	ostream& operator<< (ostream& os, const Problem& pbm)
	{
		os << endl << endl << "Number of Variables " << pbm._dimension
		   << endl;

		//Imprimo el arreglo con los costos
		os<<"Matriz de costos: "<<endl<<endl;
		for (int i=0;i<pbm._dimension;i++){
			for(int j=0;j<pbm._dimension;j++)
				os<<pbm._matrizCostos[i][j]<<" ";
			os<<endl;
		}
		os<<endl;

		//Imprimo el arreglo con los costos
		os<<"Matriz de temporadas: "<<endl<<endl;
		for (int i=0;i<3;i++){
			for(int j=0;j<2;j++)
				os<<pbm._matrizTemporadas[i][j]<<" ";
			os<<endl;
		}
		os<<endl;

		return os;
	}

	istream& operator>> (istream& is, Problem& pbm)
	{
		char buffer[MAX_BUFFER];
		int i;

		is.getline(buffer,MAX_BUFFER,'\n');
		sscanf(buffer,"%d",&pbm._dimension);

		//Pido memoria para almacenar la matriz de costos
		pbm._matrizCostos = new int *[pbm._dimension];
		for (int i=0;i<pbm._dimension;i++){
			pbm._matrizCostos[i] = new int [pbm._dimension];
			for (int j=0;j<pbm._dimension;j++)
				pbm._matrizCostos[i][j]=-1;
		}
		
            
            
            //Obtengo los nombres de las matrices asociadas
            FILE* stream = fopen("nombreMatrices.txt", "r");
            char costos[255];
            char temporadas[255];
            
            fgets(costos, 1024, stream);
            fgets(temporadas, 1024, stream);
            
            cout << "costos: " << strtok(costos,"\n") << endl;
            cout << "temporadas: " << strtok(temporadas,"\n") << endl;
            

		//Cargo el archivo con costos
	     stream = fopen(costos, "r");

	    char line[1024];
	    int j=0;
	   	while (fgets(line, 1024, stream))
	    {

	        char* tmp = strdup(line);
	        for (int i=0; i<pbm._dimension; i++){
		        const int costo = pbm.getField(tmp,i+1,' ');
				tmp = strdup(line);
		        pbm._matrizCostos[j][i] = costo;
		    }
		    j++;
		    free(tmp);
	    }

	    //Pido memoria para almacenar las temporadas
		for (int i=0;i<3;i++){
			for (int j=0;j<2;j++)
				pbm._matrizTemporadas[i][j]=0;
		}

		//Cargo el archivo con temporadas
	    stream = fopen(temporadas, "r");

	    j=0;
	    while (fgets(line, 1024, stream))
	    {

	        char* tmp = strdup(line);
	        for (int i=0; i<2; i++){
		        const int temporada = pbm.getField(tmp,i+1,',');
		        pbm._matrizTemporadas[j][i] = temporada;
		    }
		    j++;
		    free(tmp);
	    }


		cout<<pbm;
		return is;
	}

	//Funcion para leer el archivo TXT
	int Problem::getField(string input, int num, char separator){
		int iter = 1;

		for(int i=0; i < input.length();i++){
			if(input[i] == separator){
				if(num == 1){
					return atoi(input.substr(0,i).c_str());
				} else {
					iter++;
					if (iter == num){
						int j = iter + 1;
						while(input[j] != separator)
						j++;

						return atoi(input.substr(i+1, j - (i + 1)).c_str());
					}
				}
			}
		}
	}

	int ** Problem::matrizCostos() const
	{
		return _matrizCostos;
	}

	int Problem::getInicioTempBaja() const
	{
		return _matrizTemporadas[0][0];
	}
	
	int Problem::getFinTempBaja() const
	{
		return _matrizTemporadas[0][1];
	}
	
	int Problem::getInicioTempMedia() const
	{
		return _matrizTemporadas[1][0];
	}
	
	int Problem::getFinTempMedia() const
	{
		return _matrizTemporadas[1][1];
	}
	
	int Problem::getInicioTempAlta() const
	{
		return _matrizTemporadas[2][0];
	}
	
	int Problem::getFinTempAlta() const
	{
		return _matrizTemporadas[2][1];
	}
	
	bool Problem::operator== (const Problem& pbm) const
	{
		if (_dimension!=pbm.dimension()) return false;
		return true;
	}

	bool Problem::operator!= (const Problem& pbm) const
	{
		return !(*this == pbm);
	}

	Direction Problem::direction() const
	{
		//return maximize;
		return minimize;
	}

	int Problem::dimension() const
	{
		return _dimension;
	}

	Problem::~Problem()
	{
	}

	// Solution --------------------------------------------------------------

	Solution::Solution (const Problem& pbm):_pbm(pbm),_var(pbm.dimension())
	{}

	const Problem& Solution::pbm() const
	{
		return _pbm;
	}

	Solution::Solution(const Solution& sol):_pbm(sol.pbm())
	{
		*this=sol;
	}

	istream& operator>> (istream& is, Solution& sol)
	{
		for (int i=0;i<sol.pbm().dimension();i++)
			is >> sol._var[i];
		return is;
	}

	ostream& operator<< (ostream& os, const Solution& sol)
	{
		for (int i=0;i<sol.pbm().dimension();i++)
			os << " " << sol._var[i];
		return os;
	}

	NetStream& operator << (NetStream& ns, const Solution& sol)
	{
		for (int i=0;i<sol._var.size();i++)
			ns << sol._var[i];
		return ns;
	}

	NetStream& operator >> (NetStream& ns, Solution& sol)
	{
		for (int i=0;i<sol._var.size();i++)
			ns >> sol._var[i];
		return ns;
	}

 	Solution& Solution::operator= (const Solution &sol)
	{
		_var=sol._var;
		return *this;
	}

	bool Solution::operator== (const Solution& sol) const
	{
		if (sol.pbm() != _pbm) return false;
		for(int i = 0; i < _var.size(); i++)
			if(_var[i] != sol._var[i]) return false;
		return true;
	}

	bool Solution::operator!= (const Solution& sol) const
	{
		return !(*this == sol);
	}

	void Solution::initialize()
	{
		int max = _pbm.dimension();
		int aux, ind1, ind2;

		_var[0] = 0; //arranca en mvd

		//inicializo el arreglo como el camino ordenado
		for (int i=1;i<_pbm.dimension();i++)
			_var[i]= i;

		// move values to create random array
		for(int i=0;i< (max*5) ; i++)
		{
			ind1 = rand_int(1,max-1);
			ind2 = rand_int(1,max-1);

			aux = _var[ind1];
			_var[ind1] = _var[ind2];
			_var[ind2] = aux;
		}
	}

	double Solution::fitness ()
	{
		//En caso que no sea montevideo la primer ciudad
		if (_var[0] != 0) 
			return float(INT_MAX);
		
		
        double fitness = 0.0d;
        int dia = 1;
        float sobrecosto_temp_media = 0.1;
        float sobrecosto_temp_alta = 0.3;
        int** costos = _pbm.matrizCostos();
        int costo = 0;
		
		for(int i = 1; i< _pbm.dimension();i++){
			costo = costos[_var[i-1]][_var[i]];
			
			if (costo == -1)
				return float(INT_MAX);
			
			if(dia < _pbm.getFinTempBaja()){
				fitness += costo;
			}
			else if(dia < _pbm.getFinTempMedia()){
				fitness += costo + round(costo * sobrecosto_temp_media);
			}
			else{
				fitness += costo + round(costo * sobrecosto_temp_alta);
			}
			
			dia += 5;
		}

		return fitness;
	}

	char *Solution::to_String() const
	{
		return (char *)_var.get_first();
	}

	void Solution::to_Solution(char *_string_)
	{
		int *ptr=(int *)_string_;
		for (int i=0;i<_pbm.dimension();i++)
		{
			_var[i]=*ptr;
			ptr++;
		}
	}

	unsigned int Solution::size() const
	{
		return (_pbm.dimension() * sizeof(int));
	}
	
	unsigned int Solution::dimension() const
	{
		return (_pbm.dimension());
	}


	int& Solution::var(const int index)
	{
		return _var[index];
	}


	Rarray<int>& Solution::array_var()
	{
		return _var;
	}

	Solution::~Solution()
	{}

	// UserStatistics -------------------------------------------------------

	UserStatistics::UserStatistics ()
	{}

	ostream& operator<< (ostream& os, const UserStatistics& userstat)
	{
		os << "\n---------------------------------------------------------------" << endl;
		os << "                   STATISTICS OF TRIALS                   	 " << endl;
		os << "------------------------------------------------------------------" << endl;

		for (int i=0;i< userstat.result_trials.size();i++)
		{
			os << endl
			   << userstat.result_trials[i].trial
			   << "\t" << userstat.result_trials[i].best_cost_trial
			   << "\t\t" << userstat.result_trials[i].worst_cost_trial
			   << "\t\t\t" << userstat.result_trials[i].nb_evaluation_best_found_trial
			   << "\t\t\t" << userstat.result_trials[i].nb_iteration_best_found_trial
			   << "\t\t\t" << userstat.result_trials[i].time_best_found_trial
			   << "\t\t" << userstat.result_trials[i].time_spent_trial;
		}
		os << endl << "------------------------------------------------------------------" << endl;
		return os;
	}

	UserStatistics& UserStatistics::operator= (const UserStatistics& userstats)
	{
		result_trials=userstats.result_trials;
		return (*this);
	}

	void UserStatistics::update(const Solver& solver)
	{
		if( (solver.pid()!=0) || (solver.end_trial()!=true)
		  || ((solver.current_iteration()!=solver.setup().nb_evolution_steps())
		       && !terminateQ(solver.pbm(),solver,solver.setup())))
			return;

		struct user_stat *new_stat;

		if ((new_stat=(struct user_stat *)malloc(sizeof(struct user_stat)))==NULL)
			show_message(7);
		new_stat->trial         		 		 = solver.current_trial();
		new_stat->nb_evaluation_best_found_trial = solver.evaluations_best_found_in_trial();
		new_stat->nb_iteration_best_found_trial  = solver.iteration_best_found_in_trial();
		new_stat->worst_cost_trial     		 	 = solver.worst_cost_trial();
		new_stat->best_cost_trial     		 	 = solver.best_cost_trial();
		new_stat->time_best_found_trial		 	 = solver.time_best_found_trial();
		new_stat->time_spent_trial 		 		 = solver.time_spent_trial();

		result_trials.append(*new_stat);
	}

	void UserStatistics::clear()
	{
		result_trials.remove();
	}

	UserStatistics::~UserStatistics()
	{
		result_trials.remove();
	}

// Intra_operator  --------------------------------------------------------------

	Intra_Operator::Intra_Operator(const unsigned int _number_op):_number_operator(_number_op),probability(NULL)
	{}

	unsigned int Intra_Operator::number_operator() const
	{
		return _number_operator;
	}

	Intra_Operator *Intra_Operator::create(const unsigned int _number_op)
	{
		switch (_number_op)
		{
			case 0: return new Crossover;break;
			case 1: return new Mutation();break;
		}
	}

	ostream& operator<< (ostream& os, const Intra_Operator& intra)
	{
		switch (intra.number_operator())
		{
			case 0: os << (Crossover&)intra;break;
			case 1: os << (Mutation&)intra;break;
		}
		return os;
	}

	Intra_Operator::~Intra_Operator()
	{}

//  Crossover:Intra_operator -------------------------------------------------------------

	Crossover::Crossover():Intra_Operator(0)
	{
		probability = new float[1];
	}

	void Crossover::cross(Solution& sol1,Solution& sol2) const // dadas dos soluciones de la poblacion, las cruza
	{
		int i,j;
		int dim = sol1.dimension();
		
		// Copy old solutions
		Rarray<int> aux1(dim);
		aux1=sol1.array_var();
		Rarray<int> aux2(dim);
		aux2=sol2.array_var();
		int limit2=rand_int(1,dim-1);
		int limit1=rand_int(0,limit2-1);

		for (i = limit1; i < limit2; i++)
		{
			//~ cout << "i=" << i << " AUX1=" << aux1[i] << " AUX2=" << aux2[i] << endl;
			sol2.var(i) = aux1[i];
			sol1.var(i) = aux2[i];
			
		}

		for (i = 0; i < limit1; i++)
		{
			sol1.var(i) = newValue(aux1[i],limit1,limit2,aux1,aux2);
			sol2.var(i) = newValue(aux2[i],limit1,limit2,aux2,aux1);
		}

		for (i = limit2; i < dim; i++)
		{
			sol1.var(i) = newValue(aux1[i],limit1,limit2,aux1,aux2);
		 	sol2.var(i) = newValue(aux2[i],limit1,limit2,aux2,aux1);
		}

		complete(sol1);
		complete(sol2);
	}
	
	// Auxiliar function for PMX
	int Crossover::newValue(const int oldValue,const int l1,const int l2, Rarray<int> & s1, Rarray<int> & s2) const
	{
		bool fin = false;
		int nv = oldValue;
		int n = s1.size();
		bool *examinado = new bool[n];

		for(int i = 0; i < n; i++) examinado[i] = false;

		while (!fin)
		{
			fin = true;
			for(int i = l1; i < l2; i++)
				if(nv == s2[i])
				{
					if(!examinado[i])
					{
        					nv = s1[i];
						examinado[i] = true;
        					fin = false;
					}
					else	nv = -1;
        				break;
				}
		}

		delete [] examinado;

		return nv;
	}

	void Crossover::complete(Solution& s) const
	{
		int num = 0;
		int n = s.dimension();//s.size();
		int j,k;
		bool *escogido = new bool[n];

		for(int i = 0; i < n; i++) escogido[i] = false;

		for(int i = 0; i < n; i++)
		{
			if(s.var(i) != -1) escogido[s.var(i)/*-1*/] = true;
			else num++;
		}

		j =  k = 0;
		for(int i = 0; i < num; i++)
		{
			while((j < n) && (s.var(j) != -1)) j++;
			while((k < n) && escogido[k]) k++;

			s.var(j) = k+1;
		}

		delete [] escogido;
	}

	void Crossover::execute(Rarray<Solution*>& sols) const
	{
		for (int i=0;i+1<sols.size();i=i+2)
		 	if (rand01()<=probability[0]) cross(*sols[i],*sols[i+1]);
	}

	ostream& operator<< (ostream& os, const Crossover&  cross)
	{
		 os << "Crossover." << " Probability: "
                    << cross.probability[0]
		    << endl;
		 return os;
	}

	void Crossover::RefreshState(const StateCenter& _sc) const
	{
		_sc.set_contents_state_variable("_crossover_probability",(char *)probability,1,sizeof(float));
	}

	void Crossover::UpdateFromState(const StateCenter& _sc)
	{
		 unsigned long nbytes,length;
		 _sc.get_contents_state_variable("_crossover_probability",(char *)probability,nbytes,length);
	}

	void Crossover::setup(char line[MAX_BUFFER])
	{
		int op;
		sscanf(line," %d %f ",&op,&probability[0]);
		assert(probability[0]>=0);
	}

	Crossover::~Crossover()
	{
		delete [] probability;
	}

	//  Mutation: Sub_operator -------------------------------------------------------------

	Mutation::Mutation():Intra_Operator(1)
	{
		probability = new float[2];
	}

	void Mutation::mutate(Solution& sol) const
	{
		const int tam = sol.pbm().dimension();
		int swap1 = 0;
		int swap2 = 0;
		swap1 = rand_int(1,tam-1);
		swap2 = rand_int(1,tam-1);
		
		while(swap1 == swap2){
			swap1 = rand_int(1,tam-1);
			swap2 =  rand_int(1,tam-1);
		}
		
		int aux = sol.var(swap1);
		sol.var(swap1) = sol.var(swap2);
		sol.var(swap2) = aux;
	}

	void Mutation::execute(Rarray<Solution*>& sols) const
	{
		for (int i=0;i<sols.size();i++)
			if(rand01() <= probability[0])	mutate(*sols[i]);
	}

	ostream& operator<< (ostream& os, const Mutation&  mutation)
	{
		os << "Mutation." << " Probability: " << mutation.probability[0]
		   << " Probability1: " << mutation.probability[1]
		   << endl;
		return os;
	}

	void Mutation::setup(char line[MAX_BUFFER])
	{
		int op;
		sscanf(line," %d %f %f ",&op,&probability[0],&probability[1]);
		assert(probability[0]>=0);
		assert(probability[1]>=0);
	}

	void Mutation::RefreshState(const StateCenter& _sc) const
	{
		_sc.set_contents_state_variable("_mutation_probability",(char *)probability,2,sizeof(probability));
	}

	void Mutation::UpdateFromState(const StateCenter& _sc)
	{
		unsigned long nbytes,length;
		_sc.get_contents_state_variable("_mutation_probability",(char *)probability,nbytes,length);
	}

	Mutation::~Mutation()
	{
		delete [] probability;
	}

// StopCondition_1 -------------------------------------------------------------------------------------

	StopCondition_1::StopCondition_1():StopCondition()
	{}

	bool StopCondition_1::EvaluateCondition(const Problem& pbm,const Solver& solver,const SetUpParams& setup)
	{
		//return ((int)solver.best_cost_trial() == pbm.dimension());
/*		if(solver.best_cost_trial() == solver.global_best_cost()){
			//Escribo el resultado en el archivo de salida
			FILE * outputFile;

			outputFile = fopen ("solucion.out","w");
			for (int i=0;i<pbm.dimension();i++){
				fprintf (outputFile, "%d ",solver.best_solution_trial().var(i));
			}
			fclose (outputFile);

			return true;
		}*/
		return false;
		//int sol = solver

	}

	StopCondition_1::~StopCondition_1()
	{}

	//------------------------------------------------------------------------
	// Specific methods ------------------------------------------------------
	//------------------------------------------------------------------------

	bool terminateQ (const Problem& pbm, const Solver& solver,
			 const SetUpParams& setup)
	{
		StopCondition_1 stop;
		return stop.EvaluateCondition(pbm,solver,setup);
	}
}
#endif

