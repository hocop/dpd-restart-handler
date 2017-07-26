Проект для обработки restart файлов.

Можно загружать и сохранять форматы:
-- mol
-- формат Алексея Гаврилова dat (и conf)
-- уменьшенный формат sdat (весит в ~10 раз меньше за счет более оптимальной записи)

Данные загружаются в структуру Restart

Текущие возможности

Операции с Restart:
-- выделить частицы только определенного типа
-- слить две системы в одну. (задаются координаты, углы Эйлера и масштабирование на случай разных плотностей). Это нужно, чтобы собирать системы из отдельных молекул.

Анализ:
-- получить список всех цепей, где цепь - это массив номеров частиц
-- получить список всех узлов (когда к одной частице прикреплено более двух других)

Анализ ЖК:
-- директор
-- нематический параметр порядка
-- смектический параметр порядка
