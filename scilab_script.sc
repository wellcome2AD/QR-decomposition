sizes = [100 200 300 400 500]
for size = 1:len(sizes)
    // прочитать A из файла
    fid = mopen("matrixA_" + string(size) + ".txt", "r")
    if fid == -1 then
        disp("Не удалось открыть файл")
    else
        content = mgetl(fid)
        mclose(fid)
        rows = strsplit(content, " ;")
        A = []
        for j = 1:size(rows)(1)-1
            nums_array = evstr(rows(j))
            A = [A; nums_array] // добавляем как новую строку
        end
    end
    disp(A)
end