// RwR~.h

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdbool.h>
#include<time.h>
#include<string.h>
#include<float.h>

typedef struct{
    float x;
    float y;
}Vector2;

typedef struct{
    Vector2 arr[50];
    int head;
}PointList;

typedef struct{
    float delay;
    float attenuate_factor;
}DG;

void getDG(DG* a, float dist, float c){
    float d_reference = 0.05;

    a->delay = dist / c;
    //a->delay = dist / c * 5500;

    a->attenuate_factor =  d_reference/dist;
}

void sortPoint(Vector2 a[], int n, float corr){
    int i, j;
    Vector2 v;

    for(i=0; i<n; i++){
        a[i].x -=corr; a[i].y -=corr;
    }

    for(i=1; i<n; i++){
        v = a[i];
        j = i-1;
        while(j >= 0 && atan2(a[j].y, a[j].x) > atan2(v.y, v.x)){
            a[j+1] = a[j];
            j--;
        }
        a[j+1] = v;
    }

    for(i=0; i<n; i++){
        a[i].x +=corr; a[i].y +=corr;
    }
}

Vector2 makeNormal(Vector2 e1, Vector2 e2){
    Vector2 normal;

    normal.x = e2.y - e1.y;
    normal.y = e1.x - e2.x;

    return normal;
}

float dot(Vector2 a, Vector2 b){
    return a.x * b.x + a.y * b.y;
}

float distance(Vector2 a, Vector2 b){
    return sqrt((b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y));
}

// LineClipping Algorithm
int isDirect(PointList hull[], Vector2 normal[], Vector2 ss, Vector2 w){
    int i;
    float tE = 0.0f, tL = 1.0f, t, numerator, denominator;
    Vector2 Pw_Ps, direct;

    // Calclating the vector from window to soundsource
    Pw_Ps.x = ss.x - w.x;
    Pw_Ps.y = ss.y - w.y;

    for(i=0; i<hull->head; i++){
        // Calclating the vector from soundsource to edges
        direct.x = hull->arr[i].x - ss.x;
        direct.y = hull->arr[i].y - ss.y;

        numerator = dot(normal[i], direct);
        denominator = dot(normal[i], Pw_Ps);

        t =  numerator / ((-1.0) * denominator);

        // when t = nan, inequality returns false.
        if(denominator > 0 && t > tE){
            tE = t;
        }
        else if(denominator < 0 && t < tL){
            tL = t;
        }
    }

    // Direct path
    if(fabsf(tE - tL) < 0.0001){
        return 1;
    }
    // Partially disturbed
    else if(tE < tL){
        return 0;
    }
    // Completely outside
    else if(tE > tL){
        return 2;
    }
    // Not expected
    else{
        return 3;
    }
}

// Returns the side of point p with respect to line
// joining points p1 and p2.
int findSide(Vector2 p1, Vector2 p2, Vector2 p)
{
    float val = (p.y - p1.y) * (p2.x - p1.x) - (p2.y - p1.y) * (p.x - p1.x);
    if (val > 0) return 1;
    if (val < 0) return -1;
    return 0;
}

// Returns a value proportional to the distance
// between the point p and the line joining the
// points p1 and p2
float lineDist(Vector2 p1, Vector2 p2, Vector2 p)
{
    return fabs ((p.y - p1.y) * (p2.x - p1.x) - (p2.y - p1.y) * (p.x - p1.x));
}

// End points of line L are p1 and p2. Side can have value
// 1 or -1 specifying each of the parts made by the line L
void quickHull(Vector2 a[], int n, Vector2 p1, Vector2 p2, int side, PointList hull[])
{
    int i;
    int ind = -1;
    float temp = 0;
    float max_dist = 0;
    bool p1flag = false, p2flag = false;

    // Finding the point with maximum distance
    // from L and also on the specified side of L.
    for (i=0; i<n; i++)
    {
        temp = lineDist(p1, p2, a[i]);
        if(findSide(p1, p2, a[i]) == side && temp > max_dist)
        {
            ind = i;
            max_dist = temp;
        }
    }

    // If no point is found, add the end points
    // of L to the convex hull.
    // Do not allow duplication.
    if (ind == -1)
    {
        for(i=0;i<hull->head;i++){
            if(hull->arr[i].x == p1.x && hull->arr[i].y == p1.y) p1flag = true;
            if(hull->arr[i].x == p2.x && hull->arr[i].y == p2.y) p2flag = true;
            if(p1flag == true && p2flag == true) break;
        }
        if(p1flag == false){
            hull->arr[hull->head] = p1;
            hull->head++;
        }
        if(p2flag == false){
            hull->arr[hull->head] = p2;
            hull->head++;
        }
        return;
    }

    // If the point is found, recursive for the two parts divided by a[ind]
    quickHull(a, n, a[ind], p1, -findSide(a[ind], p1, p2), hull);
    quickHull(a, n, a[ind], p2, -findSide(a[ind], p2, p1), hull);
}

void buildHull(Vector2 a[], int n, PointList hull[])
{
    // Finding the point with minimum and
    // maximum x-coordinate
    int min_x = 0, max_x = 0;
    for (int i=1; i<n; i++)
    {
        if (a[i].x < a[min_x].x) min_x = i;
        if (a[i].x > a[max_x].x) max_x = i;
    }

    // Recursively find convex hull points on
    // one side of line joining a[min_x] and
    // a[max_x]
    quickHull(a, n, a[min_x], a[max_x], 1, hull);

    // Recursively find convex hull points on
    // other side of line joining a[min_x] and
    // a[max_x]
    quickHull(a, n, a[min_x], a[max_x], -1, hull);
}
