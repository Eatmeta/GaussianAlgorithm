using System;
using System.Collections.Generic;
using System.Linq;

namespace GaussAlgorithm
{
    public class Solver
    {
        private const double ZeroEdge = 1e-6;
        public double[] Solve(double[][] matrix, double[] freeMembers)
        {
            var result = new double?[matrix[0].Length];

            // Получим лист листов.
            var listList = GetListOfLists(matrix, freeMembers);

            // Если все члены уравнений равны нулю.
            if (listList.All(list => list.All(d => d == 0)))
                return new double[matrix[0].Length];

            // Если система уравнений переопредленная - удалим дублирующиеся уравнения.
            while (DeleteRedundantEquation(listList))
            {
            }

            // Решений нет если имеются равные правые части, но при этом не равные левые.
            var repeatedRightElements = listList
                .Select(list => list[list.Count - 1])
                .GroupBy(v => v).Where(g => g.Count() > 1)
                .Select(g => g.Key)
                .ToList();
            foreach (var repeatedElement in repeatedRightElements)
            {
                var leftParts = listList
                    .Where(list => Math.Abs(list[list.Count - 1] - repeatedElement) < 1e-6)
                    .Select(l => l.Take(listList[0].Count - 1).ToList())
                    .ToList();

                if (DeleteRedundantEquation(leftParts))
                    throw new NoSolutionException("System has no solution!");
            }

            // Решений нет если имеются равные левые части, но при этом не равные правые.
            if (DeleteRedundantEquation(listList.Select(list => list.Take(list.Count - 1).ToList()).ToList()))
            {
                if (listList
                        .Select(list => list[list.Count - 1])
                        .GroupBy(v => v).Where(g => g.Count() > 1)
                        .Select(g => g.Key)
                        .ToList().Count < 2)
                    throw new NoSolutionException("System has no solution!");
            }

            // Если переменных больше, чем уравнений в системе.
            var difference = listList[0].Count - listList.Count;
            if (difference > 0 && listList.Count > 1)
            {
                for (var i = listList[0].Count - 2; i > listList[0].Count - difference - 1; i--)
                {
                    result[i] = 0;
                    foreach (var list in listList)
                        list[i] = 0;
                }
            }

            // Если мы имеем систему с только двумя уравнения и только двумя ненулевыми элементами, то решим ее.
            if (listList.Count == 2 && listList.All(l => l.Take(l.Count - 1).Count(d => d != 0) == 2))
            {
                SolveSystemOfTwoEquationsWithoutZeroes(listList, result);
                return result.Select(d => (double)d).ToArray();
            }

            // Если один из элементов во всех уравнениях равен 0, то сразу пишем ему в результат 0.
            DeleteAlwaysZeroElement(listList, result);

            while (result.Any(d => d == null) || listList.Count != 0)
            {
                if (SolveSingleEquation(listList, result))
                    continue;

                // Домножаем уравнения чтобы по методу Гаусса получить нули в данном уравнении.
                MultiplyByGaussianMethod(listList);

                // Система не имеет решения, если все коэффициенты равны нулю, кроме последнего.
                if (CheckForNoSolution(listList))
                    throw new NoSolutionException("System has no solution!");

                // Вычислим корни из преобразованной системы уравнений.
                // Найдем уравнения где только один элемент в левой части не равен нулю и вычислим его значение.
                while (RemoveEquationWithOnlyOneNotZeroElement(listList, result))
                {
                }
                
                // Если образовались уравнения со сплошными нулями - удалим их.
                while (GetIndexOfEquationWithOnlyZeroes(listList) != null)
                    listList.RemoveAt((int)GetIndexOfEquationWithOnlyZeroes(listList));
            }
            return result.Select(d => (double)d).ToArray();
        }

        private bool CheckForNoSolution(List<List<double>> listList)
        {
            return listList.Any(list => list.Take(list.Count - 1).All(d => d == 0) && list.Last() != 0);
        }

        private static void MultiplyByGaussianMethod(List<List<double>> listList)
        {
            for (var i = 0; i < listList.Count; i++)
            {
                for (var j = 0; j < listList.Count; j++)
                {
                    if (i == j || listList[i][i] == 0)
                        continue;
                    var multiplier = listList[j][i] / listList[i][i];
                    // Если последний элемент превращается в ноль, то лучше вообще не трогать это уравнение.
                    if (Math.Abs(listList[j].Last() - listList[i].Last() * multiplier) < ZeroEdge)
                        continue;

                    for (var l = 0; l < listList[j].Count; l++)
                    {
                        listList[j][l] -= listList[i][l] * multiplier;
                        if (Math.Abs(listList[j][l]) < ZeroEdge && Math.Abs(listList[j][l]) > 0)
                            listList[j][l] = 0;
                    }
                }
                // Проверим, если появилось уравнение только с одним элементом не равным нулю, то выходим.
                if (listList.Any(l => l.Take(l.Count - 1).Count(e => e != 0) == 1))
                    return;
            }
        }

        private bool SolveSingleEquation(List<List<double>> listList, double?[] result)
        {
            if (listList.Count != 1) return false;
            // Если в уравнении есть больше двух элементов не равных нулю 
            if (listList[0].Take(listList[0].Count - 1).Count(e => e != 0) > 2)
            {
                var indexes = GetIndexesOfNotZeroElements(listList, 0);
                result[indexes[0]] = listList[0].Last() / listList[0][indexes[0]];
                for (var i = 1; i < indexes.Count; i++)
                    result[i] = 0;
            }
            // Если в уравнении есть два элемента не равных нулю 
            if (listList[0].Take(listList[0].Count - 1).Count(e => e != 0) == 2)
            {
                var indexes = GetIndexesOfNotZeroElements(listList, 0);
                result[indexes[0]] = 0;
                result[indexes[1]] = listList[0].Last() / listList[0][indexes[1]];
                listList.RemoveAt(0);
                return true;
            }
            // Если в уравнении нет двух элементов не равных нулю 
            if (listList[0].Take(listList[0].Count - 1).Count(e => e != 0) < 2)
            {
                for (var i = 0; i < listList[0].Count - 1; i++)
                {
                    if (listList[0][i] == 0)
                        result[i] = 0;
                    else
                        result[i] = listList[0].Last() / listList[0][i];
                }
                listList.RemoveAt(0);
                return true;
            }
            return false;
        }

        private static void DeleteAlwaysZeroElement(List<List<double>> listList, double?[] result)
        {
            for (var i = 0; i < listList[0].Count - 1; i++)
            {
                if (listList.All(list => list[i] == 0))
                    result[i] = 0;
            }
        }

        private void SolveSystemOfTwoEquationsWithoutZeroes(List<List<double>> listList, double?[] result)
        {
            var a = GetIndexesOfNotZeroElements(listList, 0);
            var b = GetIndexesOfNotZeroElements(listList, 1);
            if (!a.SequenceEqual(b)) return;

            var second = (listList[1].Last() * listList[0][a[0]] - listList[1][a[0]] * listList[0].Last()) /
                         (listList[1][a[1]] * listList[0][a[0]] - listList[1][a[0]] * listList[0][a[1]]);
            var first = (listList[0].Last() - listList[0][a[1]] * second) / listList[0][a[0]];
            result[a[0]] = first;
            result[a[1]] = second;
        }

        private List<int> GetIndexesOfNotZeroElements(List<List<double>> listList, int i)
        {
            var indexes = new List<int>();
            for (var j = 0; j < listList[i].Count - 1; j++)
            {
                if (listList[i][j] != 0)
                    indexes.Add(j);
            }
            return indexes;
        }

        private static bool RemoveEquationWithOnlyOneNotZeroElement(List<List<double>> listList, double?[] result)
        {
            foreach (var list in listList.Where(l => l.Take(l.Count - 1).Count(e => e != 0) == 1))
            {
                for (var i = 0; i < list.Count; i++)
                {
                    if (list[i] == 0) continue;
                    result[i] = list[list.Count - 1] / list[i];
                    listList.Remove(list);
                    ReplaceCalculatedElementInOtherEquations(i, (double) result[i], listList);
                    return true;
                }
            }
            return false;
        }

        private static void ReplaceCalculatedElementInOtherEquations(int i, double d, List<List<double>> listList)
        {
            foreach (var list in listList)
            {
                list[list.Count - 1] -= list[i] * d;
                list[i] = 0;
            }
        }

        private int? GetIndexOfEquationWithOnlyZeroes(List<List<double>> listList)
        {
            for (var i = 0; i < listList.Count; i++)
            {
                if (listList[i].All(d => d == 0))
                    return i;
            }
            return null;
        }

        private bool DeleteRedundantEquation(List<List<double>> listList)
        {
            var counter = 1;
            foreach (var current in listList)
            {
                foreach (var list in listList.Skip(counter))
                {
                    var sumList = current.Zip(list, (x, y) => x / y).ToList();
                    if (sumList.All(e => Math.Abs(e - sumList[0]) < ZeroEdge))
                    {
                        listList.Remove(list);
                        return true;
                    }
                }
                counter++;
            }
            return false;
        }

        private static List<List<double>> GetListOfLists(double[][] matrix, double[] freeMembers)
        {
            var listList = new List<List<double>>();
            for (var i = 0; i < freeMembers.Length; i++)
            {
                listList.Add(new List<double>());

                for (var j = 0; j < matrix[0].Length; j++)
                    listList[i].Add(matrix[i][j]);

                listList[i].Add(freeMembers[i]);
            }
            return listList;
        }
    }
}