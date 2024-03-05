#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

#include "functions.h"

#define CMAX 10000

/*** 訪問チェック関数 最近棒法用　***/
int is_in_path(int p[], int s, int n)
{
  int t, r = 0;
  for (t = 0; t < n; t++)
    /* 作成する */
    if (s == p[t])
    {
      r = 1;
    }
  return r;
}

/*** 近傍法 ***/
double solve_neighborhoods(int p[], double d[N][N], int n)
{
  double sum = 0.0, min;
  int i, j, k;
  p[0] = 0;
  for (i = 0; i < n - 1; i++)
  {
    min = d[0][0];
    for (j = 1; j < n; j++)
    {
      if (!is_in_path(p, j, i))
        if (min > d[p[i]][j])
        {
          min = d[p[i]][j];
          k = j;
        }
    }
    sum += d[p[i]][k];
    p[i + 1] = k;
  }
  sum += d[p[n - 1]][0];
  return sum;
}

void Usage(char *s)
{
  printf("Usage: %s file or num\n", s);
  exit(0);
}

// 2つの都市間のユークリッド距離を計算する関数
double distance(Location city1, Location city2)
{
  return sqrt((city1.x - city2.x) * (city1.x - city2.x) + (city1.y - city2.y) * (city1.y - city2.y));
}

// 2opt法
void solve_two_opt(int p[], double d[N][N], Location city[], int n)
{
  for (int o = 0; o < twoOPT_EPOCH; o++)
  {
    int i = 0;
    int j = 0;

    while (abs(i - j) <= 1)
    {
      i = (int)(drand48() * n);
      j = (int)(drand48() * n);
    }

    if (i > j)
    {
      int num = i;
      i = j;
      j = num;
    }

    if ((d[p[i]][p[j]] + d[p[i + 1]][p[j + 1]]) < (d[p[i]][p[i + 1]] + d[p[j]][p[j + 1]]))
    {
      for (int k = 0; k < (int)(abs(j - i) / 2); k++)
      {
        int num = p[i + 1 + k];
        p[i + 1 + k] = p[j - k];
        p[j - k] = num;
      }
    }
  }
}

/*** 経路順を表示する関数 ***/
double two_opt_show_path(int p[], double d[N][N], int n)
{
  int t;
  double sum = 0.0;
  for (t = 0; t < n; t++)
  {
    // printf("%d ", p[t]);
    sum += d[p[t]][p[t + 1]];
  }
  printf("2optのpath length = %f\n", sum);
  return sum;
}

double calc_path(int p[], double d[N][N], int n)
{
  int t;
  double sum = 0.0;
  for (t = 0; t < n; t++)
  {
    // printf("%d ", p[t]);
    sum += d[p[t]][p[t + 1]];
  }
  return sum;
}

// 経路の適応度（総距離）を計算する関数
double calc_fitness(Route population, double d[N][N], int n)
{
  int t;
  double sum = 0.0;
  for (t = 0; t < n; t++)
  {
    sum += d[population.city[t]][population.city[t + 1]];
  }
  return sum;
}

// ランダムな経路を初期化する関数
void initializeRoute(Route *city, int n)
{
  int availableCities[n]; // 利用可能な都市のリスト

  // 利用可能な都市リストを2からnまでの整数で初期化
  for (int i = 0; i < n - 1; i++)
  {
    availableCities[i] = i + 1;
  }

  city->city[0] = 0;
  city->city[n] = 0;

  for (int i = 1; i < n; i++)
  {
    // 利用可能な都市からランダムに都市を選び、経路に追加
    int randomIndex = rand() % (n - i);
    city->city[i] = availableCities[randomIndex];

    // 選択した都市を利用可能な都市のリストから削除
    for (int j = randomIndex; j < n - 1; j++)
    {
      availableCities[j] = availableCities[j + 1];
    }
  }
}

// ランキング選択２
Route ranking_selection2(Route population[POPULATION_SIZE])
{
  int selectNumber1 = rand() % POPULATION_SIZE;
  int selectNumber2 = rand() % POPULATION_SIZE;
  if (population[selectNumber1].fitness < population[selectNumber2].fitness)
  {
    return population[selectNumber1];
  }
  else
  {
    return population[selectNumber2];
  }
}

// ランキング選択3
Route ranking_selection3(Route population[POPULATION_SIZE])
{
  int selectNumber1 = rand() % POPULATION_SIZE;
  int selectNumber2 = rand() % POPULATION_SIZE;
  int selectNumber3 = rand() % POPULATION_SIZE;
  double smallest = population[selectNumber1].fitness;
  int selectNumber = selectNumber1;

  if (population[selectNumber2].fitness < smallest)
  {
    selectNumber = selectNumber2;
  }

  if (population[selectNumber3].fitness < smallest)
  {
    selectNumber = selectNumber3;
  }
  return population[selectNumber];
}

// 部分的交叉
void partical_crossover(Route *parent1, Route *parent2, int n,
                        Route *child1, Route *child2, int children_amount)
{
  double random_rate = (double)rand() / RAND_MAX;
  // printf("random_rate: %lf\n", random_rate);
  if (random_rate < CROSSOVER_RATE)
  {                                        // 交叉率によって交叉するかしないか
    int cycleStart = 1 + rand() % (n - 1); // 親1のcityの配列の番号.開始位置をランダムに選択
    int parent1_city = parent1->city[cycleStart];
    int parent2_city = parent2->city[cycleStart];
    int serch_parent1_city;
    int serch_parent2_city;
    int cross_location1 = 0;
    int cross_location2 = 0;
    do
    {
      cross_location2++;
      serch_parent2_city = parent2->city[cross_location2];
    } while (parent1_city != serch_parent2_city);

    do
    {
      cross_location1++;
      serch_parent1_city = parent1->city[cross_location1];
    } while (parent2_city != serch_parent1_city);

    child1->city[cycleStart] = parent2_city;
    child1->city[cross_location1] = parent1_city;

    child2->city[cycleStart] = parent1_city;
    child2->city[cross_location2] = parent2_city;

    for (int i = 0; i < n; i++)
    { // 親の遺伝子を交叉させるところ以外コピー
      if (i != cycleStart)
      {
        if (i != cross_location1)
        {
          child1->city[i] = parent1->city[i];
        }
        if (i != cross_location2)
        {
          child2->city[i] = parent2->city[i];
        }
      }
    }
  }
  else
  { // 交叉しないとき
    for (int i = 0; i < n; i++)
    { // 親の遺伝子をコピー
      child1->city[i] = parent1->city[i];
      child2->city[i] = parent2->city[i];
    }
  }
}

// 配列をスワップするユーティリティ関数
void swap(int *a, int *b)
{
  int temp = *a;
  *a = *b;
  *b = temp;
}

// クイックソートのパーティション関数
int partition(int array[], int low, int high)
{
  int pivot = array[high]; // ピボットを選択
  int i = low - 1;

  for (int j = low; j < high; j++)
  {
    if (array[j] < pivot)
    {
      i++;
      swap(&array[i], &array[j]);
    }
  }

  swap(&array[i + 1], &array[high]);
  return i + 1;
}

// クイックソート関数
void quickSort(int array[], int low, int high)
{
  if (low < high)
  {
    int pi = partition(array, low, high);

    quickSort(array, low, pi - 1);
    quickSort(array, pi + 1, high);
  }
}

void OX2(Route parent1, Route parent2, Route *child, int n, int children_amount)
{
  double random_rate = (double)rand() / RAND_MAX;
  // printf("random_rate: %lf\n", random_rate);
  if (random_rate < CROSSOVER_RATE)
  {
    int availableCitieslocation[n]; // 利用可能な配列の番号のリスト。
    int cycleStart[OX2location];    // ランダムな位置.配列の番号。配列の0は変更しないので１からの値。
    // 利用可能な配列の番号を1からn-1までの整数で初期化
    for (int i = 0; i < n - 1; i++)
    {
      availableCitieslocation[i] = i + 1;
    }

    for (int i = 0; i < OX2location; i++)
    {
      // 利用可能な都市からランダムに都市を選び、経路に追加
      int randomIndex = rand() % (n - i - 1);
      cycleStart[i] = availableCitieslocation[randomIndex];

      // 選択した都市を利用可能な都市のリストから削除
      for (int j = randomIndex; j < n - 1; j++)
      {
        availableCitieslocation[j] = availableCitieslocation[j + 1];
      }
    }
    // 上記は1~都市数-1までかぶらないランダムな数字をOX2locationだけ作成

    int cross_location2[OX2location]; // 　親２の場所
    for (int i = 0; i < OX2location; i++)
    { // 初期化
      cross_location2[i] = 0;
    }

    for (int i = 0; i < OX2location; i++)
    {
      do
      {
        cross_location2[i]++;
      } while (parent1.city[cycleStart[i]] != parent2.city[cross_location2[i]]);
    }

    quickSort(cross_location2, 0, OX2location - 1);
    quickSort(cycleStart, 0, OX2location - 1);

    int a = 0;
    for (int index = 0; index < n; index++)
    {
      if (index == cross_location2[a])
      {
        // 親１の要素をコピー
        child->city[index] = parent1.city[cycleStart[a]];
        a++;
      }
      else // 交叉する場所でない場合
      {
        // 親2の要素をコピー
        child->city[index] = parent2.city[index];
      }
    }
  }
  else
  { // 交叉しないとき
    for (int i = 0; i < n; i++)
    { // 親の遺伝子をコピー
      child->city[i] = parent1.city[i];
    }
  }
}

void OX(int n, Route parent1, Route parent2, Route *child)
{

  int start; // 開始点と終了点
             // ランダムな開始点を生成（1からnまで）
  start = 1 + (rand() % (n - 1));

  int divide = n / 3;
  // printf("start: %d\n", start);
  // printf("end: %d\n", end);
  // printf("n: %d\n", n);
  int crosslocationamount = 0;
  // printf("parent1の経路 ");
  // for (int i = 0; i < n; i++)
  // {
  //   printf("%d ", parent1.city[i]);
  // }
  // printf("\n");
  // printf("parent2の経路 ");
  // for (int i = 0; i < n; i++)
  // {
  //   printf("%d ", parent2.city[i]);
  // }
  // printf("\n");
  int cross_location2[N]; // 　親２の場所
  for (int i = 0; i < n; i++)
  { // 初期化
    cross_location2[i] = 0;
  }

  // printf("初めの分割");
  for (int i = 0; i < start - 1; i++)
  {
    // printf("%d", parent1.city[i + 1]);
    do
    {
      cross_location2[i]++;
    } while (parent1.city[i + 1] != parent2.city[cross_location2[i]]);
    crosslocationamount++;
  }
  // printf("\n");
  // printf("cross_location2: ");
  // for (int i = 0; i < start - 1; i++)
  // {
  //   printf("%d ", cross_location2[i]);
  // }
  // printf("\n");
  // quickSort(cross_location2, 0, start - 2); // これいらなくね？

  // printf("cross_location2: ");
  // for (int i = 0; i < start - 1; i++)
  // {
  //   printf("%d ", cross_location2[i]);
  // }
  // printf("\n");

  // printf("次の分割");
  // for (int i = start; i < start + divide; i++)
  // {
  //   printf("%d", parent1.city[i]);
  // }
  // printf("\n");

  // printf("最後の分割");
  for (int i = start + divide; i < n; i++)
  {
    // printf("%d", parent1.city[i]);
    do
    {
      cross_location2[i - divide - 1]++;
    } while (parent1.city[i] != parent2.city[cross_location2[i - divide - 1]]);
    crosslocationamount++;
  }
  // printf("\n");

  // printf("cross_location2_end: ");
  // for (int i = 0; i < crosslocationamount; i++)
  // {
  //   printf("%d ", cross_location2[i]);
  // }
  // printf("\n");
  quickSort(cross_location2, 0, crosslocationamount - 1);

  // printf("cross_location2_end: ");
  // for (int i = 0; i < crosslocationamount; i++)
  // {
  //   printf("%d ", cross_location2[i]);
  // }
  // printf("\n");

  // printf("crosslocationamount: %d\n", crosslocationamount);
  int a = 0;
  child->city[0] = 0;
  child->city[n] = 0;
  for (int index = 1; index < start; index++)
  {
    child->city[index] = parent2.city[cross_location2[a]];
    a++;
  }
  for (int index = start; index < start + divide; index++)
  {
    child->city[index] = parent1.city[index];
  }
  for (int index = start + divide; index < n; index++)
  {
    child->city[index] = parent2.city[cross_location2[a]];
    a++;
  }
  // printf("childの経路 ");
  // for (int i = 0; i <= n; i++)
  // {
  //   printf("%d ", child->city[i]);
  // }
  // printf("\n");
}

// 逆位の突然変異
void inversion_mutation(int n, Route parent, Route *child)
{
  double random_rate = (double)rand() / RAND_MAX;
  if (random_rate < MUTATION_RATE)
  {
    for (int i = 0; i < n; i++)
    { // 親の遺伝子をコピー
      child->city[i] = parent.city[i];
    }

    int start, end; // 開始点と終了点
                    // ランダムな開始点を生成（1からnまで）
    start = 1 + (rand() % (n - 1));
    end = start + (rand() % (n - start));

    while (start < end)
    {
      int temp = child->city[start];
      child->city[start] = child->city[end];
      child->city[end] = temp;
      start++;
      end--;
    }
  }
  else
  { // 交叉しないとき
    for (int i = 0; i < n; i++)
    { // 親の遺伝子をコピー
      child->city[i] = parent.city[i];
    }
  }
}

void SA_inversion_mutation(int n, Route parent, Route *child)
{
  for (int i = 0; i < n; i++)
  { // 親の遺伝子をコピー
    child->city[i] = parent.city[i];
  }

  int start, end; // 開始点と終了点
                  // ランダムな開始点を生成（1からnまで）
  start = 1 + (rand() % (n - 1));
  end = start + (rand() % (n - start));

  while (start < end)
  {
    int temp = child->city[start];
    child->city[start] = child->city[end];
    child->city[end] = temp;
    start++;
    end--;
  }
}

// 位置移動の突然変異
void locationMove_mutation(int n, Route parent, Route *child)
{
  double random_rate = (double)rand() / RAND_MAX;
  if (random_rate < MUTATION_RATE)
  {
    int front, back; // 開始点と終了点
                     // ランダムな開始点を生成（1からn-1まで）
    front = 1 + (rand() % (n - 2));
    back = front + 1 + (rand() % (n - front - 1));
    // printf("start: %d\n", front);
    // printf("end: %d\n", back);
    // printf("n: %d\n", n);
    // printf("parentの経路 ");
    // for (int i = 0; i < n; i++)
    // {
    //   printf("%d ", parent.city[i]);
    // }
    // printf("\n");
    for (int i = 0; i < n; i++)
    {
      if (i < front || i > back)
      {
        child->city[i] = parent.city[i];
      }
      else if (i < back)
      {
        child->city[i + 1] = parent.city[i];
      }
    }

    child->city[front] = parent.city[back];
    child->city[front + 1] = parent.city[front];
    // printf("childの経路 ");
    // for (int i = 0; i <= n; i++)
    // {
    //   printf("%d ", child->city[i]);
    // }
    // printf("\n");
  }
  else
  { // 交叉しないとき
    for (int i = 0; i < n; i++)
    { // 親の遺伝子をコピー
      child->city[i] = parent.city[i];
    }
  }
}

Route maxORmin_fitness(Route population[POPULATION_SIZE], bool max)
{
  Route elite = population[0];
  for (int i = 0; i < POPULATION_SIZE; i++)
  {
    if (max)
    {
      if (population[i].fitness < elite.fitness)
      {
        elite = population[i];
      }
    }
    else
    {
      if (population[i].fitness > elite.fitness)
      {
        elite = population[i];
      }
    }
  }
  return elite;
}

double avg_fitness(Route population[POPULATION_SIZE])
{
  double sum = 0.0;
  for (int i = 0; i < POPULATION_SIZE; i++)
  {
    sum += population[i].fitness;
  }
  return sum / POPULATION_SIZE;
}

// 配列を比較し、同じ場合は1、異なる場合は0を返す関数
int areArraysEqual(int arr1[], int arr2[], int size)
{
  for (int i = 0; i < size; i++)
  {
    if (arr1[i] != arr2[i])
    {
      return 0; // 配列が異なる
    }
  }
  return 1; // 配列が同じ
}

// 配列を表示する関数
void printArray(int arr[], int size)
{
  for (int i = 0; i < size; i++)
  {
    printf("%d ", arr[i]);
  }
  printf("\n");
}

// 焼きなまし法
void simulatedAnnealing(int initial[], int n, double d[N][N])
{
  time_t startTime, currentTime;
  double elapsedTime;
  int completedSA = 0;
  // プログラム開始時刻を記録
  time(&startTime);

  Route bestPath[SATIMES];
  Route elite;

  for (int k = 0; k < SATIMES; k++)
  {
    Route currentPath;

    Route newCurrentPath;
    double temperature = INITIAL_TEMPERATURE;

    // Copy the initial path to current path and best path
    for (int i = 0; i < n + 1; i++)
    {
      currentPath.city[i] = bestPath[k].city[i] = initial[i];
    }
    // printf("currentPath\n");
    // printArray(currentPath.city,n+1);

    bestPath[k].fitness = currentPath.fitness = calc_path(initial, d, n);

    while (temperature > FINAL_TEMPERATURE)
    {
      SA_inversion_mutation(n, currentPath, &newCurrentPath);

      double newDistance = calc_fitness(newCurrentPath, d, n);

      double delta = newDistance - currentPath.fitness;

      // printf("currentPath: %f, newDistance: %f\n", currentPath.fitness, newDistance);

      if (delta <= 0 || rand() / (double)RAND_MAX <= exp(-delta / temperature))
      {
        currentPath.fitness = newDistance;
        if (newDistance < bestPath[k].fitness)
        {
          bestPath[k].fitness = newDistance;
          for (int i = 0; i < n + 1; i++)
          {
            bestPath[k].city[i] = newCurrentPath.city[i];
            currentPath.city[i] = newCurrentPath.city[i];
          }
        }
      }

      temperature *= COOLING_RATE;
    }
    completedSA++;
    if (completedSA % elapsedTimeNumberSA == 0)
    {
      elite = bestPath[0];
      for (int i = 0; i < completedSA; i++)
      {
        if (bestPath[i].fitness < elite.fitness)
        {
          elite = bestPath[i];
        }
      }
       printf("n: %d\n", n);
      int n_acc = 0;
      for (int f = 0; f < n; f++)
      {
        n_acc += is_in_path(elite.city, f, n);
      }
      printf("n_acc: %d\n", n_acc);
      if(n != n_acc) {
        exit(0);
      }
      printf("%d 回の一番良い適応度: %lf\n", completedSA, elite.fitness);
      time(&currentTime);
      elapsedTime = difftime(currentTime, startTime);
      double remainingTime = (elapsedTime / completedSA) * (SATIMES - completedSA);

      printf("進行中... 現在の経過時間: %.2f分, 残り時間: %.2f分\n", elapsedTime / 60, remainingTime / 60);
    }
  }

  elite = bestPath[0];
  for (int i = 0; i < SATIMES; i++)
  {
    if (bestPath[i].fitness < elite.fitness)
    {
      elite = bestPath[i];
    }
  }
  for (int i = 0; i < n + 1; i++)
  {
    printf("%d ", elite.city[i]);
    initial[i] = elite.city[i];
  }

  printf("最終的な焼きなましのfitness: %f\n", elite.fitness);
  time(&currentTime);
  elapsedTime = difftime(currentTime, startTime);
  printf("SA処理完了。総経過時間: %.2f分\n", elapsedTime / 60);
}

// クイックソートのパーティション関数
int Routepartition(Route array[], int low, int high)
{
  double pivot = array[high].fitness;
  int i = low - 1;

  for (int j = low; j < high; j++)
  {
    if (array[j].fitness < pivot)
    {
      i++;
      // Route 構造体の要素を交換
      Route temp = array[i];
      array[i] = array[j];
      array[j] = temp;
    }
  }

  // Route 構造体の要素を交換
  Route temp = array[i + 1];
  array[i + 1] = array[high];
  array[high] = temp;

  return i + 1;
}

// クイックソート関数
void RoutequickSort(Route array[], int low, int high)
{
  if (low < high)
  {
    int pi = Routepartition(array, low, high);
    RoutequickSort(array, low, pi - 1);
    RoutequickSort(array, pi + 1, high);
  }
}

int main(int argc, char *argv[])
{
  int n, i, j;
  Location location[N];
  double dist[N][N], len;
  int path[N];
  bool last = false;               // 最終結果を表示するため
  bool init = true;                // fitnessfileの初期化用
  srand((unsigned int)time(NULL)); // rand用

  char input_path_file[50] = pathFile;
  char input_fitness_file[50] = fitnessFile;

  if (argc < 2)
    Usage(argv[0]);
  /** ノード数 n **/
  n = atoi(argv[1]);
  if (n > 0)
  {
    n = atoi(argv[1]);
    // random_location(n, location);
  }
  else
  {
    n = load_location(argv[1], location);
    if (CROSS == "partical")
    { // 部分交叉の場合ファイル名にpartを付ける
      strcat(input_path_file, "part_");
      strcat(input_fitness_file, "part_");
    }
    else if (CROSS == "OX2")
    {
      strcat(input_path_file, "OX2_");
      strcat(input_fitness_file, "OX2_");
    }
    else if (CROSS == "OX")
    {
      strcat(input_path_file, "OX_");
      strcat(input_fitness_file, "OX_");
    }
    else
    {
      printf("綴りが間違ってるよ\n");
    }
    // ファイル名に入力したファイル名を追加
    strcat(input_path_file, argv[1]);
    strcat(input_fitness_file, argv[1]);
    save_city_fitness(input_fitness_file, 0.0, last, init);
    init = false;
  }

  /*** 距離行列の計算 ***/
  calc_distance(location, dist, n);

  /*** 近傍法で解く ***/
  len = solve_neighborhoods(path, dist, n);
  printf("最近傍法の path length = %f\n", len);
  save_city_fitness(input_fitness_file, len, last, init);

  // 焼きなまし法で解く
  if (isSA)
  {
    printf("SA開始\n");
    simulatedAnnealing(path, n, dist);
  }
  printf("n: %d\n", n);
  int n_acc = 0;
  for (int f = 0; f < n; f++)
  {
    n_acc += is_in_path(path, f, n);
  }
  printf("n_acc: %d\n", n_acc);
  if(n != n_acc) {
        exit(0);
      }

  time_t startTime, currentTime;
  double elapsedTime;
  int completedGenerations = 0;

  // プログラム開始時刻を記録
  time(&startTime);

  printf("GA開始\n");
  printf("ELITE: %d\n", ELITE);
  printf("SELECT_CHILD: %d\n", SELECT_CHILD);
  printf("PARENT_ALIVE_RATE: %f\n", PARENT_ALIVE_RATE);
  printf("CROSSOVER_RATE: %f\n", CROSSOVER_RATE);
  printf("MUTATION_RATE: %f\n", MUTATION_RATE);
  printf("GENERATION: %d\n", GENERATION);
  printf("CROSS: %s\n", CROSS);
  printf("SELECTION: %d\n", RANKING_SELECTION);
  printf("MUTATION: %s\n", MUTATION);
  printf("people_twoOPT: %d\n", people_twoOPT);

  // 初期候補経路を生成
  Route population[POPULATION_SIZE];
  Route child[SELECT_CHILD];
  Route elite;

  // 　初期個体の生成
  for (int i = 0; i < POPULATION_SIZE; i++)
  {
    // printf("%d 個目の初期経路の都市名の順番: ",i);
    initializeRoute(&population[i], n);
  }

  readNumbersFromFile(bestscoreFile, path, N);

  for (int i = 0; i < ELITE; i++)
  {
    for (int j = 0; j < n; j++)
    {
      population[POPULATION_SIZE - 1 - i].city[j] = path[j];
    }
  }

  // 世代ループ
  for (int i = 0; i < GENERATION; i++)
  {
    // 各個体の評価
    for (int j = 0; j < POPULATION_SIZE; j++)
    {
      population[j].fitness = calc_fitness(population[j], dist, n);
    }

    // エリート選択のために一番適応度が低い個体を抽出
    elite = maxORmin_fitness(population, true);
    Route dropouter = maxORmin_fitness(population, false);
    double avgfitness = avg_fitness(population);

    int increase_child = 1;
    if (CROSS == "partical")
    {
      increase_child = 2;
    }

    // 子の生成
    for (int children_amount = (int)(SELECT_CHILD * PARENT_ALIVE_RATE); children_amount < SELECT_CHILD;
         children_amount += increase_child) // particalの場合子が2体ずつ、ox2の場合1体ずつ
    {
      Route parent1;
      Route parent2;
      if (RANKING_SELECTION == 2)
      {
        // ランキング選択2
        parent1 = ranking_selection2(population);
        parent2 = ranking_selection2(population);
      }
      else
      {
        // ランキング選択3
        parent1 = ranking_selection3(population);
        parent2 = ranking_selection3(population);
      }

      // 同じ親をとらないようにする
      while (areArraysEqual(parent1.city, parent2.city, n))
      {
        // 一方の配列を変更
        parent2 = ranking_selection2(population);
      }

      if (CROSS == "partical")
      {
        // 　部分的交叉
        partical_crossover(&parent1, &parent2, n, &child[children_amount], &child[children_amount + 1], children_amount);
        inversion_mutation(n, child[children_amount + 1], &child[children_amount + 1]); // 逆位の突然変異
      }
      else if (CROSS == "OX2")
      {
        // 一様順序交叉
        OX2(parent1, parent2, &child[children_amount], n, children_amount);
      }
      else if (CROSS == "OX")
      {
        OX(n, parent1, parent2, &child[children_amount]);
      }
      else
      {
        printf("綴りが間違ってるよ\n");
      }

      if (MUTATION == "inversion")
      {
        // 逆位の突然変異.親を突然変異させる
        inversion_mutation(n, child[children_amount], &child[children_amount]);
      }
      else if (MUTATION == "locationMove")
      {
        locationMove_mutation(n, child[children_amount], &child[children_amount]);
      }
      else if (MUTATION == "SA")
      {
        // simulatedAnnealing(child[children_amount].city,n,dist);
      }
      else
      {
        printf("綴りが間違ってるよ\n");
      }
    }

    RoutequickSort(population, 0, POPULATION_SIZE - 1);

    // 世代交代
    for (int j = (int)(SELECT_CHILD * PARENT_ALIVE_RATE); j < SELECT_CHILD; j++)
    {
      population[j] = child[j];
    }

    // 最後の子にエリートを設定
    for (int i = 0; i < ELITE; i++)
    {
      population[POPULATION_SIZE - i - ELITE] = elite;
    }

    for (int l = 0; l < people_twoOPT; l++)
    {
      int randam = (rand() % (POPULATION_SIZE - 1));
      solve_two_opt(population[randam].city, dist, location, n);
    }

    completedGenerations++; // 世代カウント

    if (completedGenerations % elapsedTimeNumber == 0)
    {
      printf("n: %d\n", n);
      int n_acc = 0;
      for (int f = 0; f < n; f++)
      {
        n_acc += is_in_path(elite.city, f, n);
      }
      printf("n_acc: %d\n", n_acc);
      printf("%d 世代のエリートの適応度: %lf, 落第者の適応度: %lf, 平均の適応度: %lf\n", completedGenerations, elite.fitness, dropouter.fitness, avgfitness);
      if(n != n_acc) {
        exit(0);
      }
        save_path_GA(input_path_file, elite.city, n, last);
        save_city_fitness(input_fitness_file, elite.fitness, last, init);
      // elapsedTimeNumber回ごとに経過時間を計算して表示
      time(&currentTime);
      elapsedTime = difftime(currentTime, startTime);
      double remainingTime = (elapsedTime / completedGenerations) * (GENERATION - completedGenerations);

      printf("進行中... 現在の経過時間: %.2f分, 残り時間: %.2f分\n", elapsedTime / 60, remainingTime / 60);
    }
  }

  printf("n: %d\n", n);
  n_acc = 0;
  for (int f = 0; f < n; f++)
  {
    n_acc += is_in_path(elite.city, f, n);
  }
  printf("n_acc: %d\n", n_acc);
  if(n != n_acc) {
        exit(0);
      }

  last = true;
  // 最後に求めたエリートを保存
  elite.fitness = calc_fitness(elite, dist, n);

  if (n == n_acc)
  {
    // 経路を保存するファイルの名前を変更
    char elite_fitness[20];
    sprintf(elite_fitness, "%.4f", elite.fitness);
    strcat(input_path_file, elite_fitness);
    save_path_GA(input_path_file, elite.city, n, last);
    strcat(input_fitness_file, elite_fitness);
    save_city_fitness(input_fitness_file, elite.fitness, last, init);
  }
  // char new_input_fitness_file[50];
  // strcpy(new_input_fitness_file, fitnessFile);
  // // strcat(new_input_fitness_file, input_fitness_file);
  // strcat(new_input_fitness_file, elite_fitness);
  // // const char *sourceFileName = input_fitness_file;
  // // const char *destinationFileName = new_input_fitness_file;
  // // printf("%s,%s\n",sourceFileName,destinationFileName);
  // copyFile(input_fitness_file, new_input_fitness_file);

  time(&currentTime);
  elapsedTime = difftime(currentTime, startTime);
  printf("GA処理完了。総経過時間: %.2f分\n", elapsedTime / 60);

  return 0;
}