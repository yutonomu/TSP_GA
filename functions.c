#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "functions.h"

/***** ノードの座標を読み込む *****/
int load_location(char *fname, Location loc[])
{
  FILE *fp;
  int n = 0, i;
  double x, y;
  if ((fp = fopen(fname, "r")) == NULL)
  {
    printf("Unable to open file %s\n", fname);
    exit(0);
  }

  while (fscanf(fp, "%d%lf%lf", &i, &x, &y) == 3)
  {
    loc[n].x = x;
    loc[n].y = y;
    n++;
  }
  fclose(fp);
  return n;
}

/***** パスの保存 *****/
void save_path(char *fname, int p[], int n)
{
  FILE *fp;
  int t;
  if ((fp = fopen(fname, "w")) == NULL)
  {
    printf("Unable to open file %s\n", fname);
    exit(0);
  }
  for (t = 0; t < n; t++)
    fprintf(fp, "%d ", p[t]);
  fprintf(fp, "%d\n", p[0]);
  fclose(fp);
}
/***** ノードの座標を出力する *****/
void save_location(char *fname, Location loc[], int p[], int n)
{
  FILE *fp;
  int t;
  if ((fp = fopen(fname, "w")) == NULL)
  {
    printf("Unable to open file %s\n", fname);
    exit(0);
  }
  for (t = 0; t < n; t++)
    fprintf(fp, "%f %f\n", loc[p[t]].x, loc[p[t]].y);
  fprintf(fp, "%f %f\n", loc[p[0]].x, loc[p[0]].y);
  fclose(fp);
}

/***** ノード間の距離を計算する *****/
void calc_distance(Location loc[], double d[N][N], int n)
{
  int i, j;
  double dist, max = 0;
  for (i = 0; i < n; i++)
  {
    for (j = i + 1; j < n; j++)
    {
      dist = sqrt((loc[i].x - loc[j].x) * (loc[i].x - loc[j].x) + (loc[i].y - loc[j].y) * (loc[i].y - loc[j].y));
      d[i][j] = dist;
      d[j][i] = dist;
      if (dist > max)
        max = dist;
    }
  }
  for (i = 0; i < n; i++)
    d[i][i] = max + 1;
  // printf("# max= %f\n", max);
}

/***** GAパスの保存 *****/
void save_path_GA(char *fname, int p[], int n, bool last)
{
  // 相対パスでの書き込み先ディレクトリ
  const char *GApathfileDir = "GApathfile/";
  char filePath[256];
  FILE *fp;
  // ファイルパスを構築
  snprintf(filePath, sizeof(filePath), "%s%s", GApathfileDir, fname);

  int t;
  if ((fp = fopen(filePath, "w")) == NULL)
  {
    printf("Unable to open file %s\n", fname);
    exit(0);
  }
  for (t = 0; t < n; t++)
    fprintf(fp, "%d ", p[t]);
  fprintf(fp, "%d\n", p[0]);
  if (last)
  {
    printf("最終的なエリートの経路: ");
    for (t = 0; t < n; t++)
      printf("%d ", p[t]);
    printf("%d\n", p[0]);
  }

  fclose(fp);
}

/***** GAの適応度の保存 *****/
void save_city_fitness(char *fname, double p, bool last, bool init)
{
  const char *fitnessfileDir = "fitnessfile/";
  char filePath[256];
  FILE *fp;
  // ファイルパスを構築
  snprintf(filePath, sizeof(filePath), "%s%s", fitnessfileDir, fname);
  if (init)
  {
    // ファイルを書き込みモードでオープンし、既存の内容を削除
    if ((fp = fopen(filePath, "w")) == NULL)
    {
      perror("ファイルをオープンできません");
      return;
    }
    // ファイルをすぐに閉じる
    fclose(fp);
  }
  else
  {
    if ((fp = fopen(filePath, "a")) == NULL)
    {
      printf("Unable to open file %s\n", fname);
      exit(0);
    }
    fprintf(fp, "%f\n", p);
    if (last)
    {
      printf("最終的なエリートの適応度: %f\n", p);
    }
    fclose(fp);
  }
}

void copyFile(const char *sourceFileName, const char *destinationFileName)
{
  const char *fitnessfileDir = "fitnessfile/";
  char sourcefilePath[256];
  char destinationfilePath[256];
  FILE *sourceFile, *destinationFile;
  // ファイルパスを構築
  snprintf(sourcefilePath, sizeof(sourcefilePath), "%s%s", fitnessfileDir, sourceFileName);
  snprintf(destinationfilePath, sizeof(destinationfilePath), "%s%s", fitnessfileDir, destinationFileName);
  char ch;

  // 既存のファイルを読み込む
  sourceFile = fopen(sourcefilePath, "rb");
  if (sourceFile == NULL)
  {
    perror("既存のファイルを開けません");
    return;
  }

  // 新しいファイルを作成または上書きモードで開く
  destinationFile = fopen(destinationfilePath, "wb");
  if (destinationFile == NULL)
  {
    perror("新しいファイルを開けません");
    fclose(sourceFile); // ソースファイルを閉じる
    return;
  }

  // 既存のファイルから新しいファイルに内容をコピー
  while ((ch = fgetc(sourceFile) != EOF))
  {
    fputc(ch, destinationFile);
    printf("ch: %d\n",ch);
  }

  // ファイルを閉じる
  fclose(sourceFile);
  fclose(destinationFile);
}

// ファイルから数字を読み込んで配列に格納する関数
int readNumbersFromFile(const char *filename, int path[], int maxPathLength) {
   const char *GApassfileDir = "GApathfile/";
  char filePath[256];
  // ファイルパスを構築
  snprintf(filePath, sizeof(filePath), "%s%s", GApassfileDir, filename);
    FILE *file = fopen(filePath, "r");
  
  
    if (file == NULL) {
        perror("ファイルを開けません");
        return 0;
    }

    int pathLength = 0;

    while (pathLength < maxPathLength && fscanf(file, "%d", &path[pathLength]) == 1) {
        pathLength++;
    }

    fclose(file);
    return pathLength;
}