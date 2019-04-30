#include <initializer_list>
#include <iostream>
#include <vector>

struct relation {
   int len;
   int* terms;
   relation(): len(0), terms(nullptr) {}
   relation(int* terms, int size) {
      initialize(terms, size);
   }
   relation(std::initializer_list<int> il_terms) {
      int* terms = new int[il_terms.size()];
      int i = 0;
      for (auto it = std::begin(il_terms); it != std::end(il_terms); ++i, ++it) {
         terms[i] = *it;
      }
      initialize(terms, il_terms.size());
      delete[] terms;
   }
   private:

   void initialize(int* terms, int size) {
      int mult = 1;
      int last = terms[size-1];
      if (last < 0) {
         mult = -last;
         size--;
      }
      this->len = size*mult;
      this->terms = new int[this->len];
      int i = 0;
      for (int j = 0; j < mult; j++) {
         for (int k = 0; k < size; k++, i++) {
            this->terms[i] = terms[k];
         }
      }
   }
};

template<int N>
struct group {
   static const int ngens = N;
   int nrels;
   relation* rels;
   group(): nrels(0), rels(nullptr) {}
   group(std::initializer_list<relation> rels): nrels(rels.size()) {
      this->rels = new relation[this->nrels];
      int i = 0;
      for (auto it = std::begin(rels); it != std::end(rels); ++it) {
         this->rels[i++] = *it;
      }
   }
   ~group() {
      delete[] rels;
   }
};

template<int N>
struct coxeter : group<N> {
   coxeter(): group<N>() {}
   coxeter(relation* rels, int size) {
      initialize(rels, size);
   }
   coxeter(std::initializer_list<relation> il_rels) {
      relation* rels = new relation[il_rels.size()];
      int i = 0;
      for (auto it = std::begin(il_rels); it != std::end(il_rels); ++i, ++it) {
         rels[i] = *it;
      }
      initialize(rels, il_rels.size());
      delete[] rels;
   }
   private:

   void initialize(relation* rels, int size) {
      this->nrels = ( N * (N - 1) ) / 2;
      if (this->nrels == 0) {
         this->rels = nullptr;
         return;
      }
      this->rels = new relation[this->nrels];
      bool* hasrel[N];
      for (int i = 0; i < N; i++) {
         hasrel[i] = new bool[i+1];
         for (int j = 0; j <= i; j++) {
            hasrel[i][j] = false;
         }
      }
      int p = 0;
      relation r;
      for (int i = 0; i < size; i++) {
         r = rels[i];
         int f = r.terms[0];
         int t = r.terms[1];
         if (t == f)
            continue;
         this->rels[p++] = r;
         if (f < t)
            std::swap(f,t);
         hasrel[f][t] = true;
      }
      for (int i = 0; i < N; i++) {
         for (int j = 0; j <= i; j++) {
            if (i != j and !hasrel[i][j])
               this->rels[p++] = {j, i, -2};
         }
      }
      for (int i = 0; i < N; i++) {
         delete[] hasrel[i];
      }
   }
};

template<int N>
struct cosets {
   std::vector<int> gens[N];
   int num_cosets;
   int next_ent_cs;
   int next_ent_gn;
   int alloc_size;
   cosets(): num_cosets(1), next_ent_cs(0), next_ent_gn(0), alloc_size(32) {
      for (int i=0; i<N; i++) {
         this->gens[i] = std::vector<int>(this->alloc_size,-1);
      }
   }
   ~cosets() {
      delete[] this->gens;
   }
   int apply(int coset, int gen) {
      return this->gens[gen][coset];
   }
   void set(int coset, int gen, int target) {
      this->gens[gen][coset] = target;
      this->gens[gen][target] = coset;
   }
   void find_next_entry() {
      while (this->next_ent_cs < this->num_cosets and this->apply(this->next_ent_cs, this->next_ent_gn) > -1) {
         this->next_ent_gn++;
         if (this->next_ent_gn >= N) {
            this->next_ent_gn %= N;
            this->next_ent_cs++;
         }
      }
   }
   bool add_coset() {
      this->find_next_entry();
      if (this->next_ent_cs >= this->num_cosets)
         return false;
      int next_coset = this->num_cosets;
      this->num_cosets++;
      if (this->alloc_size < this->num_cosets) {
         this->alloc_size <<= 1;
         for (int gen = 0; gen < N; gen++) {
            this->gens[gen].resize(this->alloc_size, -1);
         }
      }
      this->set(this->next_ent_cs, this->next_ent_gn, next_coset);
      return true;
   }
   void print() {
      std::cout<<"    ";
      for (int j = 0; j < N; j++) {
         std::cout<<j<<' ';
      }
      std::cout<<std::endl;
      for (int i = 0; i < this->num_cosets; i++) {
         std::cout<<i<<" | ";
         for (int j = 0; j < N; j++) {
            std::cout<<this->apply(i,j)<<' ';
         }
         std::cout<<std::endl;
      }
   }
};

template<int N>
int** cosets_to_arr(cosets<N>* c) {
   int** arr = new int*[c->num_cosets];
   for (int i = 0; i < c->num_cosets; i++) {
      int* row = new int[N];
      arr[i] = row;
      for (int j = 0; j < N; j++) {
         row[j] = c->apply(i, j);
      }
   }
   return arr;
}

struct reltable_row {
   int left_coset;
   int right_target;
   int* left_gen;
   int* right_gen;
   reltable_row* next;
   reltable_row(): left_coset(0), right_target(0), left_gen(nullptr), right_gen(nullptr), next(nullptr) {}
   reltable_row(relation* rel, int i): next(nullptr) {
      this->left_coset = i;
      this->right_target = i;
      this->left_gen = &(rel->terms[0]);
      this->right_gen = &(rel->terms[rel->len-1]);
   }
   template<int N>
   bool apply_learn(cosets<N>* c) {
      while (this->left_gen < this->right_gen) {
         int left_target = c->apply(this->left_coset, *this->left_gen);
         if (left_target < 0)
            break;
         this->left_gen++;
         this->left_coset = left_target;
      }
      while (this->left_gen < this->right_gen) {
         int right_coset = c->apply(this->right_target, *this->right_gen);
         if (right_coset < 0)
            break;
         this->right_gen--;
         this->right_target = right_coset;
      }
      if (this->right_gen == this->left_gen) {
         c->set(this->left_coset, *this->left_gen, this->right_target);
         return true;
      }
      return false;
   }
};

template<int N>
struct reltable {
   group<N>* g;
   reltable_row* first;
   reltable_row* last;
   reltable_row* unused;
   int nrows;
   reltable(group<N>* g): g(g), first(nullptr), last(nullptr), unused(nullptr), nrows(0) {}
   ~reltable() {
      reltable_row* curr, *next;
      curr = this->first;
      while (curr != nullptr) {
         next = curr->next;
         delete curr;
         curr = next;
      }
      curr = this->unused;
      while (curr != nullptr) {
         next = curr->next;
         delete curr;
         curr = next;
      }
   }
   void add_rows() {
      int rownum = 0;
      reltable_row* row;
      if (this->first == nullptr) {
         if (this->unused != nullptr) {
            row = this->unused;
            *row = reltable_row(&(this->g->rels[0]), this->nrows);
            this->unused = this->unused->next;
         }
         else {
            row = new reltable_row(&(this->g->rels[0]), this->nrows);
         }
         this->first = row;
         rownum++;
      }
      if (rownum == 0) {
         if (this->unused != nullptr) {
            row = this->unused;
            *row = reltable_row(&(this->g->rels[0]), this->nrows);
            this->unused = this->unused->next;
         }
         else {
            row = new reltable_row(&(this->g->rels[0]), this->nrows);
         }
         this->last->next = row;
         rownum++;
      }
      while (this->unused != nullptr and rownum < this->g->nrels) {
         row->next = this->unused;
         *(row->next) = reltable_row(&(this->g->rels[rownum]), this->nrows);
         this->unused = this->unused->next;
         row = row->next;
         rownum++;
      }
      for (int i = rownum; i < this->g->nrels; i++) {
         row->next = new reltable_row(&(this->g->rels[i]), this->nrows);
         row = row->next;
      }
      this->last=row;
      this->nrows++;
   }
   bool apply_learn(cosets<N>* c) {
      bool learned = false;
      reltable_row* next;
      while (this->first != nullptr and this->first->apply_learn(c)) {
         learned = true;
         next = this->first->next;
         this->first->next = this->unused;
         this->unused = this->first;
         this->first = next;
      }
      if (this->first == nullptr) {
         this->last = nullptr;
         return learned;
      }
      reltable_row* curr = this->first;
      while (curr->next != nullptr) {
         next = curr->next;
         if (next->apply_learn(c)) {
            learned = true;
            curr->next = next->next;
            next->next = this->unused;
            this->unused = next;
         }
         else
            curr = curr->next;
      }
      this->last = curr;
      return learned;
   }
};

template<int N>
cosets<N>* enumerate_cosets(group<N>* g, std::initializer_list<int> gens) {
   cosets<N>* c = new cosets<N>();
   reltable<N> rts(g);

   for (auto it = std::begin(gens); it != std::end(gens); ++it)
      c->set(0, *it, 0);

   do {
      rts.add_rows();
      while (rts.apply_learn(c));
   } while(c->add_coset());
   
   return c;
}

int main(){
   //coxeter<6> g = {{0,1,-3},{1,2,-3},{2,3,-3},{2,4,-3},{4,5,-3}};
   coxeter<7> g = {{0,1,-3},{1,2,-3},{2,3,-3},{2,4,-3},{4,5,-3},{5,6,-3}};
   //coxeter<4> g = {{0,1,-5},{1,2,-3},{2,3,-3}};
   auto c = enumerate_cosets(&g,{0,1,2});
   std::cout<<c->num_cosets<<std::endl; 
}
