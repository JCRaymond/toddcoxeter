#include <initializer_list>
#include <iostream>
#include <vector>
#include <string>

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

struct path_elem {
   int from;
   int gen;
   path_elem(): from(-1), gen(-1) {}
   path_elem(int from, int gen): from(from), gen(gen) {}
};

struct word {
   int* chars;
   int len;
   word(): chars(nullptr), len(0) {}
   word(int* chars, int len): chars(chars), len(len) {}
   word append(word other) {
      int newlen = this->len + other.len;
      int* newchars = new int[newlen];
      int p = 0;
      for (int i = 0; i < this->len; i++) {
         newchars[p++] = this->chars[i];
      }
      for (int i = 0; i < other.len; i++) {
         newchars[p++] = other.chars[i];
      }
   }
};

#define DEFAULT_ALLOC_SIZE 32

template<int N>
struct cosets {
   std::vector<path_elem> path;
   std::vector<int> gens[N];
   int num_cosets;
   int next_ent_cs;
   int next_ent_gn;
   int alloc_size;
   cosets(): num_cosets(1), next_ent_cs(0), next_ent_gn(0), alloc_size(DEFAULT_ALLOC_SIZE) {
      this->path = std::vector<path_elem>(this->alloc_size);
      this->path[0] = path_elem(0, -1);
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
   int apply(int coset, word w) {
      for (int i = 0; i < w.len; i++) {
         coset = this->apply(coset, w.chars[i]);
      }
      return coset;
   }
   int apply(int coset, word w, int* gen_map) {
      for (int i = 0; i < w.len; i++) {
         coset = this->apply(coset, gen_map[w.chars[i]]);
      }
      return coset;
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
      if (this->next_ent_cs >= this->num_cosets) {
         this->path.resize(this->num_cosets);
         for (int gen = 0; gen < N; gen++) {
            this->gens[gen].resize(this->num_cosets);
         }
         return false;
      }
      int next_coset = this->num_cosets;
      this->num_cosets++;
      if (this->alloc_size < this->num_cosets) {
         this->alloc_size <<= 1;
         for (int gen = 0; gen < N; gen++) {
            this->gens[gen].resize(this->alloc_size, -1);
         }
         this->path.resize(this->alloc_size);
      }
      this->path[next_coset] = path_elem(this->next_ent_cs, this->next_ent_gn);
      this->set(this->next_ent_cs, this->next_ent_gn, next_coset);
      return true;
   }
   word* path_to_words(int* gen_map = nullptr) {
      word* words = new word[this->num_cosets];
      for (int i = 0; i < this->num_cosets; i++) {
         words[i] = word();
      }
      int size, curr;
      for (int i = this->num_cosets - 1; i > 0; i--) {
         if (words[i].len != 0)
            continue;
         size = 0;
         curr = i;
         while (curr != 0) {
            size++;
            curr = this->path[curr].from;
         }
         int* chars = new int[size];
         size--;
         curr = i;
         if (gen_map == nullptr) {
            while (curr != 0) {
               words[curr] = word(chars, size+1);
               chars[size--] = this->path[curr].gen;
               curr = this->path[curr].from;
            }
         }
         else {
            while (curr != 0) {
               words[curr] = word(chars, size+1);
               chars[size--] = gen_map[this->path[curr].gen];
               curr = this->path[curr].from;
            } 
         }
      }
      return words;
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
      if (this->g->nrels == 0)
         return;
      int rownum = 0;
      reltable_row* row;
      if (this->unused != nullptr) {
         row = this->unused;
         this->unused = this->unused->next;
         *row = reltable_row(&(this->g->rels[0]), this->nrows);
      }
      else {
         row = new reltable_row(&(this->g->rels[0]), this->nrows);
      }
      rownum++;
      if (this->first == nullptr)
         this->first = row;
      else
         this->last->next = row;
      while (this->unused != nullptr and rownum < this->g->nrels) {
         row->next = this->unused;
         this->unused = this->unused->next;
         *(row->next) = reltable_row(&(this->g->rels[rownum]), this->nrows);
         row = row->next;
         rownum++;
      }
      for (; rownum < this->g->nrels; rownum++) {
         row->next = new reltable_row(&(this->g->rels[rownum]), this->nrows);
         row = row->next;
      }
      this->last = row;
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
struct sub_coxeter;

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
   cosets<N>* enumerate_cosets(std::initializer_list<int> gens) {
      if (N == 0)
         return new cosets<N>();
      cosets<N>* c = new cosets<N>();
      reltable<N> rts(this);
      for (auto it = std::begin(gens); it != std::end(gens); ++it) {
         c->set(0, *it, 0);
      }

      do {
         rts.add_rows();
         while (rts.apply_learn(c));
      } while(c->add_coset());
      
      return c;
   }
   cosets<N>* enumerate_cosets(int* gens, int n) {
      if (N == 0)
         return new cosets<N>();
      cosets<N>* c = new cosets<N>();
      reltable<N> rts(this);

      for (int i = 0; i < n; i++) {
         c->set(0, gens[i], 0);
      }

      do {
         rts.add_rows();
         while (rts.apply_learn(c));
      } while(c->add_coset());
      
      return c;
   }
   cosets<N>** enumerate_all_cosets() {
      int num_subgroups = 1 << N;
      int temparr[N];
      int n, mask; 
      cosets<N>** ec = new cosets<N>*[num_subgroups];
      for (int i = 0; i < num_subgroups; i++) {
         n = 0;
         mask = 0;
         for (int j = 0; j < N; j++) {
            mask = 1 << j;
            if (i & mask) {
               temparr[n++] = j;
            }
         }
         ec[i] = this->enumerate_cosets(temparr, n);
      }
      return ec;
   }
   template<int M>
   sub_coxeter<M>* get_subgroup(int* gens) {
      int nrels = ( M * (M - 1) ) / 2;
      relation* rels = new relation[nrels];
      int gen_map[N];
      for (int i = 0; i < N; i++) {
         gen_map[i] = -1;
      }
      for (int i = 0; i < M; i++) {
         gen_map[gens[i]] = i;
      }
      int p = 0;
      relation r;
      int f,t;
      for (int i = 0; i < this->nrels; i++) {
         r = this->rels[i];
         f = gen_map[r.terms[0]];
         t = gen_map[r.terms[1]];
         if (f != -1 and t != -1) {
            rels[p++] = {f, t, -r.len/2};
         }
      }
      return new sub_coxeter<M>(rels, nrels, gens);
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
struct sub_coxeter : coxeter<N> {
   int* gen_map;
   sub_coxeter(): coxeter<N>() {}
   sub_coxeter(relation* rels, int size, int* gen_map): coxeter<N>(rels, size), gen_map(gen_map) {}
};

void print_words(word* words, int n, int* gen_map = nullptr) {
   word w;
   for (int i = 0; i < n; i++) {
      w = words[i];
      std::cout<<i<<": ";
      for (int j = 0; j < w.len; j++) {
         if (gen_map != nullptr)
            std::cout<<gen_map[w.chars[j]];
         else
            std::cout<<w.chars[j];
      }
      std::cout<<std::endl;
   }
}

int main(){
   //coxeter<6> g = {{0,1,-3},{1,2,-3},{2,3,-3},{2,4,-3},{4,5,-3}};
   //coxeter<7> g = {{0,1,-3},{1,2,-3},{2,3,-3},{2,4,-3},{4,5,-3},{5,6,-3}};
   coxeter<4> g = {{0,1,-5},{1,2,-3},{2,3,-3}};
   //coxeter<3> g = {{0,1,-4},{1,2,-3}};
   int gl[] = {0,1,2};
   const int n = 2;
   auto sg = g.get_subgroup<n>(&gl[0]);

   auto t_c = g.enumerate_cosets({});
  
   print_words(t_c->path_to_words(), t_c->num_cosets);
   std::cout<<std::endl;

   auto sg_c = g.enumerate_cosets(gl, n);
   auto sg_c_words = sg_c->path_to_words();
   print_words(sg_c_words, sg_c->num_cosets);
   std::cout<<std::endl;

   auto e_sg = sg->enumerate_cosets({});
   auto e_sg_words = e_sg->path_to_words(sg->gen_map);
   print_words(e_sg_words, e_sg->num_cosets);
   std::cout<<std::endl;

   int face_points[t_c->num_cosets];
   int coset;
   for (int i = 0; i < e_sg->num_cosets; i++) {
      coset = t_c->apply(0, e_sg_words[i]);
      for (int j = 0; j < sg_c->num_cosets; j++) {
         face_points[t_c->apply(coset, sg_c_words[j])] = j;
      }
   }

   for (int i = 0; i < t_c->num_cosets; i++) {
      std::cout<<i<<": "<<face_points[i]<<std::endl;
   }
   std::cout<<std::endl;
}
