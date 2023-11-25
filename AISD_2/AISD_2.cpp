#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
using namespace std;


//Bubble sort
void BubbleSort(int arr[], int size) {
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                swap(arr[j], arr[j + 1]);
            }
        }
    }
}
//Selection sort
void SelectionSort(int a[], int N) {
    int min = 0;
    for (int i = 0; i < N; i++)
    {
        min = i;
        for (int j = i + 1; j < N; j++)
            min = (a[j] < a[min]) ? j : min;
        if (i != min)
        {
            swap(a[i], a[min]);
        }
    }
}
//Insertion sort
void InsertionSort(int arr[], int n)
{
    int i, key, j;
    for (i = 1; i < n; i++) {
        key = arr[i];
        j = i - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}
//Merge sort
void merge(int arr[], int p, int q, int r) {

    int n1 = q - p + 1;
    int n2 = r - q;

    int* L = new int[n1]; int* M = new int[n1];

    for (int i = 0; i < n1; i++)
        L[i] = arr[p + i];
    for (int j = 0; j < n2; j++)
        M[j] = arr[q + 1 + j];

    int i, j, k;
    i = 0;
    j = 0;
    k = p;

    while (i < n1 && j < n2) {
        if (L[i] <= M[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = M[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = M[j];
        j++;
        k++;
    }
}

void mergeSort(int arr[], int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;

        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        merge(arr, l, m, r);
    }
}
//Quick sort
int partition(int arr[], int start, int end)
{

    int pivot = arr[start];

    int count = 0;
    for (int i = start + 1; i <= end; i++) {
        if (arr[i] <= pivot)
            count++;
    }

    int pivotIndex = start + count;
    swap(arr[pivotIndex], arr[start]);

    int i = start, j = end;

    while (i < pivotIndex && j > pivotIndex) {

        while (arr[i] <= pivot) {
            i++;
        }

        while (arr[j] > pivot) {
            j--;
        }

        if (i < pivotIndex && j > pivotIndex) {
            swap(arr[i++], arr[j--]);
            
        }
    }

    return pivotIndex;
}

void quickSort(int arr[], int start, int end)
{

    if (start >= end)
        return;

    int p = partition(arr, start, end);
    
    quickSort(arr, start, p - 1);

    quickSort(arr, p + 1, end);

}
//Heap sort
void heapify(int arr[], int N, int i)
{

    int largest = i;

    int l = 2 * i + 1;

    int r = 2 * i + 2;

    if (l < N && arr[l] > arr[largest])
        largest = l;

    if (r < N && arr[r] > arr[largest])
        largest = r;

    if (largest != i) {
        swap(arr[i], arr[largest]);

        heapify(arr, N, largest);
    }
}

void heapSort(int arr[], int N)
{

    for (int i = N / 2 - 1; i >= 0; i--)
        heapify(arr, N, i);

    for (int i = N - 1; i > 0; i--) {

        swap(arr[0], arr[i]);

        heapify(arr, i, 0);
    }
}
//Tim sort
void insertionSort_for_tim(int arr[], int left, int right)
{
    for (int i = left + 1; i <= right; i++) {
        int temp = arr[i];
        int j = i - 1;
        while (j >= left && arr[j] > temp) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = temp;
    }
}

void timSort(int arr[], int n)
{
    int RUN = 32;
    
    for (int i = 0; i < n; i += RUN)
        insertionSort_for_tim(arr, i, min((i + RUN - 1), (n - 1)));

    
    for (int size = RUN; size < n; size = 2 * size) {

        
        for (int left = 0; left < n; left += 2 * size) {

            
            int mid = left + size - 1;
            int right = min((left + 2 * size - 1), (n - 1));

            
            if (mid < right)
                merge(arr, left, mid, right);
        }
    }
}
//Introsort
void swapValue(int* a, int* b)
{
    int* temp = a;
    a = b;
    b = temp;
    return;
}

void InsertionSort_for_introsort(int arr[], int* begin, int* end)
{
    
    int left = begin - arr;
    int right = end - arr;

    for (int i = left + 1; i <= right; i++)
    {
        int key = arr[i];
        int j = i - 1;

        
        while (j >= left && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }

    return;
}

int* Partition_for_introsort(int arr[], int low, int high)
{
    int pivot = arr[high];    
    int i = (low - 1);   

    for (int j = low; j <= high - 1; j++)
    {
        
        if (arr[j] <= pivot)
        {
             
            i++;

            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return (arr + i + 1);
}

int* MedianOfThree(int* a, int* b, int* c)
{
    if (*a < *b && *b < *c)
        return (b);

    if (*a < *c && *c <= *b)
        return (c);

    if (*b <= *a && *a < *c)
        return (a);

    if (*b < *c && *c <= *a)
        return (c);

    if (*c <= *a && *a < *b)
        return (a);

    if (*c <= *b && *b <= *a)
        return (b);
}

void IntrosortUtil(int arr[], int* begin, int* end, int depthLimit)
{
     
    int size = end - begin;

     
    if (size < 16)
    {
        InsertionSort_for_introsort(arr, begin, end);
        return;
    }

    
    if (depthLimit == 0)
    {
        make_heap(begin, end + 1);
        sort_heap(begin, end + 1);
        return;
    }

    
    int* pivot = MedianOfThree(begin, begin + size / 2, end);

     
    swapValue(pivot, end);

     
    int* partitionPoint = Partition_for_introsort(arr, begin - arr, end - arr);
    IntrosortUtil(arr, begin, partitionPoint - 1, depthLimit - 1);
    IntrosortUtil(arr, partitionPoint + 1, end, depthLimit - 1);

    return;
}

void Introsort(int arr[], int* begin, int* end)
{
    int depthLimit = 2 * log(end - begin);

     
    IntrosortUtil(arr, begin, end, depthLimit);

    return;
}
//Shell sort
int shellSort(int arr[], int n)
{
    
    for (int gap = n / 2; gap > 0; gap /= 2)
    {
        
        for (int i = gap; i < n; i += 1)
        {
            
            int temp = arr[i];

             
            int j;
            for (j = i; j >= gap && arr[j - gap] > temp; j -= gap)
                arr[j] = arr[j - gap];

             
            arr[j] = temp;
        }
    }
    return 0;
}

void ArrayIsSorted(int a[], int size) {
    int flag = 0;
    for (int i = 0; i < size - 1; i++) {
        if (a[i] > a[i + 1]) {
            cout << endl << "Array is not sorted" << endl;
            flag = 1;
            return;
        }
        else flag = 0;
    }
    cout << endl << "Array is sorted" << endl;
}

void reverce(int a[], int size) {
    for (int i = 0; i < size / 2; i++) {
        swap(a[i], a[size - 1 - i]);
    }
}

int compare(const void* a, const void* b)
{
    const int* x = (int*)a;
    const int* y = (int*)b;

    if (*x > *y)
        return 1;
    else if (*x < *y)
        return -1;

    return 0;
}

void main()
{
    fstream f;
    f.open("graph.py", ios::out);
    f << "from matplotlib import pyplot as plt" << endl;
    f << "from scipy.optimize import curve_fit\n"
        << "import numpy as np\n\n";
    f << "def f(x, a) :\n\treturn (a * x * np.log(x))\n";
    //f << "def f(x, a) :\n\treturn (a * x)\n";
    //f << "def f(x, a) :\n\treturn (a * x * x)\n";
    f << "x = np.array([";
    const int n = 50000;
    const int delta = 250;
    int k = 0;
    double* time = new double[n / delta];
    for (int i = 0; i <= n; i += delta) {
        int size = i;
        if (i == 0) size = 1;
        int *a = new int[size];
        random_device dev;
        mt19937 rng(dev());
        for (int i = 0; i < size; i++) {
            uniform_int_distribution<std::mt19937::result_type> dist6(1, 1000);
            a[i] = dist6(rng);
        }
        if (i != n) f << size << " ,";
        else if (i == n) f << size << "])";
        //timSort(a, size/2);
        //reverce(a, size);
        ArrayIsSorted(a, size); cout << endl;
        //int start, end;
        

        auto start = chrono::steady_clock::now();
        Introsort(a, a, a + size - 1);

        auto end = chrono::steady_clock::now();
        chrono::duration<double> tm = end - start;
        time[i / delta] = tm.count();
        
        
        
        
        k++;
        ArrayIsSorted(a, size); 
        cout << endl << k;
        cout << endl << (time[i / delta]);
        delete[] a;
    }
    f << endl << "y = np.array([";
    for (int i = 0; i <= n; i += delta) {
        if (i != n) f << time[i / delta] << " ,";
        else f << time[i / delta] <<  "])";
    }
    
    f << "\narg, _ = curve_fit(f, x, y)";
    f << "\nprint(arg)";
    f << "\nfunc = arg * x * np.log(x)";
    //f << "\nfunc = arg * x";
    //f << "\nfunc = arg * x * x";
    f << "\nplt.xlabel(\"Number of elements\")\nplt.ylabel(\"Time\")";
    f << endl << "plt.title(\"Semi-sorted_IntroSort\")";
    f << endl << "plt.plot(x,func, linewidth=1.5)";
    f << endl << "plt.scatter(x,y, c='red', s = 0.3)";
    //f << endl << "plt.savefig(\"./plots/Semi-sorted_IntroSort.png\")";
    f << endl << "plt.show()";
    f.close();
    system("python graph.py");
}
