// $Header$
/* ----------------------------------------------------------------------------
  selector.h
  mbwall 10aug94
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  The selectors are functions that use information in fitness objects to pick
genomes from a population.  There are a number of selection criteria, and
not every selection method will work with every fitness object.

  These are the pre-defined selection methods.  A brief outline of what each 
one does is given here, but for more details (especially the theoretical 
foundations for each method), see Goldberg.  See the source file for some
implementation details that Goldberg does not delve into.
  The selection methods can, to some extent, be mixed and matched with the
various scaling methods.
  If you want to do speciating selection (e.g. DeJong-style crowding with
Hamming distance) then you will have to write your own selection method.

         Rank - pick the genome with the best fitness (not objective score)
RouletteWheel - weighted selection where individuals with better fitness have
                a greater chance of being selected than those with lower 
                scores.
   Tournament - similar to roulette, but instead of choosing one, choose two
                then pick the better of the two as the selected individuals
      Uniform - stochastic uniform selection picks randomly from the 
                population.  Each individual has as much chance as any other.
          SRS - stochastic remainder selection does a preselection based on 
                the expected number of each genome, then a random sampling on 
                the preselected list.
           DS - deterministic sampling is implemented as described in 
                Goldberg's book (as much as I could understand it, anyway).
---------------------------------------------------------------------------- */
#ifndef _ga_selector_h_
#define _ga_selector_h_

#include <string.h>
#include <ga/gaid.h>

class GAGenome;
class GAPopulation;

/* ----------------------------------------------------------------------------
   The base class definition for the selector object defines the interface.
Any derived selector must define a clone member and a select member.  If you
add any special data members then you should also define a copy member.
   Any selector can do its business based on fitness or objective scores.  The
base selector provides the mechanism for this.  Derived classes can use it if
they want to, or ignore it.
---------------------------------------------------------------------------- */
class GASelectionScheme  : public GAID {
public:
  GADefineIdentity("GASelectionScheme", GAID::Selection);
  enum { RAW, SCALED };

  GASelectionScheme(int w=SCALED) { which = w;}
  GASelectionScheme(const GASelectionScheme& orig) { copy(orig); }
  GASelectionScheme& operator=(const GASelectionScheme& orig)
    { if(&orig != this) copy(orig); return *this; }
  virtual ~GASelectionScheme() {}
  virtual GASelectionScheme* clone() const=0;
  virtual void copy(const GASelectionScheme& orig)
    { pop=orig.pop; which=orig.which; }
  virtual void assign(GAPopulation& p) { pop = &p; }
  virtual void update() {}
  virtual GAGenome& select() const=0;

protected:
  GAPopulation* pop;
  int which;			// should we use fitness or objective scores?
};


/* ----------------------------------------------------------------------------
   The rank selector simply picks the best individual in the population.  You 
can specify whether the selector should use the raw (objective) scores or the
scaled (fitness) scores to determine the best individual.  Default is fitness.
---------------------------------------------------------------------------- */
#if USE_RANK_SELECTOR == 1
class GARankSelector : public GASelectionScheme {
public:
  GADefineIdentity("GARankSelector", GAID::RankSelection);

  GARankSelector(int w=GASelectionScheme::SCALED) : GASelectionScheme(w) {}
  GARankSelector(const GARankSelector& orig) { copy(orig); }
  GARankSelector& operator=(const GASelectionScheme& orig) 
    { if(&orig != this) copy(orig); return *this; }
  virtual ~GARankSelector() {}
  virtual GASelectionScheme* clone() const { return new GARankSelector; }
  virtual GAGenome& select() const;
};
#endif

  
/* ----------------------------------------------------------------------------
   Roulette wheel uses a fitness-proportional algorithm for selecting 
individuals.
---------------------------------------------------------------------------- */
#if USE_ROULETTE_SELECTOR == 1 || USE_TOURNAMENT_SELECTOR == 1 
class GARouletteWheelSelector : public GASelectionScheme {
public:
  GADefineIdentity("GARouletteWheelSelector", GAID::RouletteWheelSelection);

  GARouletteWheelSelector(int w=GASelectionScheme::SCALED) : 
  GASelectionScheme(w) 
    { psum = (float*)0; n = 0; }
  GARouletteWheelSelector(const GARouletteWheelSelector& orig) 
    { psum = (float*)0; n = 0; copy(orig); }
  GARouletteWheelSelector& operator=(const GASelectionScheme& orig) 
    { if(&orig != this) copy(orig); return *this; }
  virtual ~GARouletteWheelSelector() { delete [] psum; }
  virtual GASelectionScheme* clone() const
    { return new GARouletteWheelSelector; }
  virtual void copy(const GASelectionScheme& orig) {
    GASelectionScheme::copy(orig);
    const GARouletteWheelSelector& sel = 
      DYN_CAST(const GARouletteWheelSelector&,orig);
    delete [] psum;
    n = sel.n; 
    psum = new float[n];
    memcpy(psum, sel.psum, n * sizeof(float));
  }
  virtual GAGenome& select() const;
  virtual void update();

protected:
  int n;
  float* psum;
};
#endif

  
/* ----------------------------------------------------------------------------
   This version of the tournament selector does two roulette wheel selections
then picks the better of the two.  We derive from the roulette wheel class so
that we can use its update method.
---------------------------------------------------------------------------- */
#if USE_TOURNAMENT_SELECTOR == 1
class GATournamentSelector : public GARouletteWheelSelector {
public:
  GADefineIdentity("GATournamentSelector", GAID::TournamentSelection);

  GATournamentSelector(int w=GASelectionScheme::SCALED) : 
  GARouletteWheelSelector(w) {}
  GATournamentSelector(const GATournamentSelector& orig) { copy(orig); }
  GATournamentSelector& operator=(const GASelectionScheme& orig) 
    { if(&orig != this) copy(orig); return *this; }
  virtual ~GATournamentSelector() {}
  virtual GASelectionScheme* clone() const
    { return new GATournamentSelector; }
  virtual GAGenome& select() const;
};
#endif


/* ----------------------------------------------------------------------------
	This version of the selector uses Lexicase selection.
---------------------------------------------------------------------------- */
#if USE_LEXICASE_SELECTOR == 1
typedef void (*DataShuffler)();
typedef float (*LexicaseEval)(GAGenome &, int indIdx);

class GALexicaseSelector : public GASelectionScheme {
public:
  GADefineIdentity("GALexicaseSelector", GAID::LexicaseSelection);

  GALexicaseSelector(int w=GASelectionScheme::SCALED) : GASelectionScheme(w) { }
  GALexicaseSelector(DataShuffler df, LexicaseEval ef, int dsize, int w=GASelectionScheme::SCALED) : 
    GASelectionScheme(w){ dshuffler=df; lexieval=ef; datasize=dsize;}
  
  GALexicaseSelector(const GALexicaseSelector& orig) { copy(orig); }
  GALexicaseSelector& operator=(const GASelectionScheme& orig) 
    { if(&orig != this) copy(orig); return *this; }
  virtual ~GALexicaseSelector() {}
  virtual GASelectionScheme* clone() const
    { return new GALexicaseSelector(dshuffler, lexieval, datasize); }
  virtual GAGenome& select() const;
  virtual void copy(const GASelectionScheme& orig) {
  	GASelectionScheme::copy(orig);
  	const GALexicaseSelector& sel = DYN_CAST(const GALexicaseSelector&, orig);
  	dshuffler = sel.dshuffler;
  	lexieval = sel.lexieval;
  	datasize = sel.datasize;
  }
  LexicaseEval Lexievaluator() const {return lexieval;}
  LexicaseEval Lexievaluator(LexicaseEval f){return (lexieval=f);}
  DataShuffler Datashuffler() const {return dshuffler;}
  DataShuffler Datashuffler(DataShuffler f);
  void SetDataSize(int s){datasize=s;}
  int GetDataSize(){return datasize;}
  
  protected:
    DataShuffler dshuffler;
    LexicaseEval lexieval;
    int datasize;
};

#endif


/* ----------------------------------------------------------------------------
	This version of the selector uses Epsilon-Lexicase selection.
---------------------------------------------------------------------------- */
#if USE_EPSILON_LEXICASE_SELECTOR == 1
#include<vector>
typedef void (*EpsilonDataShuffler)();
typedef float (*EpsilonLexicaseEval)(GAGenome &, int indIdx);
typedef float (*EpsilonLexicaseScore)(GAGenome &, int indIdx);

class GAEpsilonLexicaseSelector : public GASelectionScheme {
public:
  GADefineIdentity("GAEpsilonLexicaseSelector", GAID::EpsilonLexicaseSelection);

  GAEpsilonLexicaseSelector(int w=GASelectionScheme::SCALED) : GASelectionScheme(w) { }
  GAEpsilonLexicaseSelector(EpsilonDataShuffler df, EpsilonLexicaseEval ef, EpsilonLexicaseScore of, int dsize, int w=GASelectionScheme::SCALED) : 
    GASelectionScheme(w){ dshuffler=df; lexieval=ef; scoreorig=of; datasize=dsize;}
  
  GAEpsilonLexicaseSelector(const GAEpsilonLexicaseSelector& orig) { copy(orig); }
  GAEpsilonLexicaseSelector& operator=(const GASelectionScheme& orig) 
    { if(&orig != this) copy(orig); return *this; }
  virtual ~GAEpsilonLexicaseSelector() {}
  virtual GASelectionScheme* clone() const
    { return new GAEpsilonLexicaseSelector(dshuffler, lexieval, scoreorig, datasize); }
  virtual GAGenome& select() const;
  virtual void copy(const GASelectionScheme& orig) {
  	GASelectionScheme::copy(orig);
  	const GAEpsilonLexicaseSelector& sel = DYN_CAST(const GAEpsilonLexicaseSelector&, orig);
  	dshuffler = sel.dshuffler;
  	lexieval = sel.lexieval;
  	scoreorig = sel.scoreorig;
  	datasize = sel.datasize;
  }
  EpsilonLexicaseEval Lexievaluator() const {return lexieval;}
  EpsilonLexicaseEval Lexievaluator(EpsilonLexicaseEval f){return (lexieval=f);}
  EpsilonLexicaseScore Epsilonlexicasescore() const {return scoreorig;}
  EpsilonLexicaseScore Epsilonlexicasescore(EpsilonLexicaseScore f){return (scoreorig=f);}
  EpsilonDataShuffler Datashuffler() const {return dshuffler;}
  EpsilonDataShuffler Datashuffler(EpsilonDataShuffler f);
  void SetDataSize(int s){datasize=s;}
  int GetDataSize(){return datasize;}
  
  protected:
  	float getMedian(std::vector<float>& ) const;
    EpsilonDataShuffler dshuffler;
    EpsilonLexicaseEval lexieval;
    EpsilonLexicaseScore scoreorig;
    int datasize;
};

#endif



/* ----------------------------------------------------------------------------
   This version of the selector uses two rounds of selection. 
 * The first round involves creating a pool of individuals selected on the 
 * basis of fitness or an alternative score in the genome.
 * The second round involves 
 * We derive from the roulette wheel class so
that we can use its update method.
---------------------------------------------------------------------------- */
#ifdef ATHENA_BLOAT_CONTROL
#if USE_DOUBLETOURNAMENT_SELECTOR == 1
class GADoubleTournamentSelector : public GASelectionScheme {
public:
  GADefineIdentity("GADoubleTournamentSelector", GAID::DoubleTournamentSelection);

  GADoubleTournamentSelector(float dval, int fval, bool ff,int w=GASelectionScheme::SCALED) : 
    GASelectionScheme(w){ D=dval; F=fval; fitness_first=ff;}
  GADoubleTournamentSelector(const GADoubleTournamentSelector& orig) { copy(orig); }
  GADoubleTournamentSelector& operator=(const GASelectionScheme& orig) 
    { if(&orig != this) copy(orig); return *this; }
  virtual ~GADoubleTournamentSelector() {}
  virtual GASelectionScheme* clone() const
    { return new GADoubleTournamentSelector(D,F,fitness_first); }
  virtual GAGenome& select() const;
  virtual int roundone_select() const;
  virtual int select_best(int ind1, int ind2, bool use_fitness) const;
  virtual void SetD(float dval){D=dval;}
  virtual void SetF(int fval){F=fval;}
  virtual void setFitnessFirst(bool ff){fitness_first=ff;}
  
  protected:
      bool fitness_first;
      float D;
      int F;
    
};

#endif
#endif


/* ----------------------------------------------------------------------------
   Stochastic uniform sampling selection.  This is just a fancy name for 
random sampling.  Any individual in the population has as much chance of being
selected as any other.
---------------------------------------------------------------------------- */
#if USE_UNIFORM_SELECTOR == 1
class GAUniformSelector : public GASelectionScheme {
public:
  GADefineIdentity("GAUniformSelector", GAID::UniformSelection);

  GAUniformSelector(int w=GASelectionScheme::SCALED) : GASelectionScheme(w) { }
  GAUniformSelector(const GAUniformSelector& orig) { copy(orig); }
  GAUniformSelector& operator=(const GASelectionScheme& orig) 
    { if(&orig != this) copy(orig); return *this; }
  virtual ~GAUniformSelector() { }
  virtual GASelectionScheme* clone() const { return new GAUniformSelector; }
  virtual GAGenome& select() const;
};
#endif


/* ----------------------------------------------------------------------------
   Stochastic remainder sampling selection.  
---------------------------------------------------------------------------- */
#if USE_SRS_SELECTOR == 1
class GASRSSelector : public GASelectionScheme {
public:
  GADefineIdentity("GASRSSelector", GAID::SRSSelection);

  GASRSSelector(int w=GASelectionScheme::SCALED) : GASelectionScheme(w)
    { fraction = (float*)0; choices = (unsigned int *)0; n = 0; }
  GASRSSelector(const GASRSSelector& orig)
    { fraction = (float*)0; choices = (unsigned int *)0; n = 0; copy(orig); }
  GASRSSelector& operator=(const GASelectionScheme& orig) 
    { if(&orig != this) copy(orig); return *this; }
  virtual ~GASRSSelector() { delete [] fraction; delete [] choices; }
  virtual GASelectionScheme* clone() const { return new GASRSSelector; }
  virtual void copy(const GASelectionScheme& orig) {
    GASelectionScheme::copy(orig);
    const GASRSSelector& sel = DYN_CAST(const GASRSSelector&, orig);
    delete [] fraction;  delete [] choices;
    n = sel.n; 
    fraction = new float [n];
    choices = new unsigned int [n];
    memcpy(fraction, sel.fraction, n * sizeof(float));
    memcpy(choices, sel.choices, n * sizeof(unsigned int));
  }
  virtual GAGenome& select() const;
  virtual void update();

protected:
  float *fraction;
  unsigned int *choices;
  unsigned int n;
};
#endif


/* ----------------------------------------------------------------------------
   Deterministic sampling selection.  
---------------------------------------------------------------------------- */
#if USE_DS_SELECTOR == 1
class GADSSelector : public GASelectionScheme {
public:
  GADefineIdentity("GADSSelector", GAID::DSSelection);

  GADSSelector(int w=GASelectionScheme::SCALED) : GASelectionScheme(w) {
    fraction = (float*)0;
    choices = (unsigned int *)0; 
    idx = (unsigned int *)0; 
    n = 0; 
  }
  GADSSelector(const GADSSelector& orig) { 
    fraction = (float*)0; 
    choices = (unsigned int *)0; 
    idx = (unsigned int *)0; 
    n = 0; 
    copy(orig); 
  }
  GADSSelector& operator=(const GASelectionScheme& orig) 
    { if(&orig != this) copy(orig); return *this; }
  virtual ~GADSSelector() {
    delete [] fraction;
    delete [] choices;
    delete [] idx;
  }
  virtual GASelectionScheme* clone() const { return new GADSSelector; }
  virtual void copy(const GASelectionScheme& orig) {
    GASelectionScheme::copy(orig);
    const GADSSelector& sel = DYN_CAST(const GADSSelector&, orig);
    delete [] fraction;  delete [] choices; delete [] idx;
    n = sel.n; 
    fraction = new float [n];
    choices = new unsigned int [n];
    idx = new unsigned int [n];
    memcpy(fraction, sel.fraction, n * sizeof(float));
    memcpy(choices, sel.choices, n * sizeof(unsigned int));
    memcpy(idx, sel.idx, n * sizeof(unsigned int));
  }
  virtual GAGenome& select() const;
  virtual void update();

protected:
  float *fraction;
  unsigned int *choices;
  unsigned int *idx;
  unsigned int n;
};
#endif


#endif
