#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#define N 480 // 都市数の上限

// 個体数を定義
#define ELITE 5 // 次の世代に残すエリートの数.
 // 交叉や突然変異（交叉では親２から子２の場合）をさせる親の数。もしくはその親から生成する子の数
#define SELECT_CHILD 1000 //　偶数にすること
#define POPULATION_SIZE (ELITE+SELECT_CHILD) //全個体数

// パラメーター
#define CROSSOVER_RATE 1.0 // 交叉率
#define MUTATION_RATE 0.3 // 突然変異率
#define PARENT_ALIVE_RATE 0.15 // 親の残存率
#define GENERATION 20000 // 世代数
#define elapsedTimeNumber 500 // 途中経過用
#define twoOPT_EPOCH 100000 // 2opt法の回数
#define people_twoOPT 6 // 2optで変化させる親の数

// SAのパラメーター
#define isSA 1 // SAを使うか。1の時使う
#define SATIMES 30
#define INITIAL_TEMPERATURE 1000000.0
#define COOLING_RATE 0.99999
#define FINAL_TEMPERATURE 0.001
#define elapsedTimeNumberSA 5 // 途中経過用

#define CROSS "OX" // どの交叉を用いるか。partical,OX2,OX
#define OX2location 4 // 一様順序交叉のランダムな位置の数.OX2の時だけ

#define RANKING_SELECTION 2 // ランキング選択2を採用の時 2

#define MUTATION "inversion" // 位置移動の突然変異の場合locationMove,逆位の時、inversion

// 以下は文字列連結によりファイル名を決定。入れるファイルの名前に気を付けること
#define pathFile "GApath_"
#define fitnessFile "fitness_"

#define bestscoreFile "bestscor" // 前のベストスコアの経路を最初のエリートに追加

typedef struct {
   double x;
   double y;
} Location;

// 候補経路を表す構造体
typedef struct {
    int city[N];
    double fitness;
} Route;

int load_location( char *fname, Location loc[] );
void load_path( char *fname, int *p, int n);
void save_location( char *fname, Location loc[], int p[], int n);
void save_path( char *fname, int p[], int n);
void save_path_GA( char *fname, int p[], int n,bool last);
void save_city_fitness( char *fname, double p,bool last,bool init);
void calc_distance( Location loc[], double d[N][N], int n );
double path_length(int *path, double d[N][N], int n);

double distance(Location city1, Location city2);
void initializeRoute(Route *city,int n);
void calculateFitness(Route *city, Location cities[N],int n);

void clearFileContents(const char *filename);

void copyFile(const char *sourceFileName, const char *destinationFileName);

int readNumbersFromFile(const char *filename, int path[], int maxPathLength);
#endif
