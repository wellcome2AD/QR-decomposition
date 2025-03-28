function matrix=readMatrixFromFile(fileName)
  fid = mopen(fileName, "r")
  if fid == -1 then
    disp("Не удалось открыть файл")
  else
    content = mgetl(fid)
    mclose(fid)
    rows = strsplit(content, " ;")
    matrix = []
    for j = 1:size(rows)(1)-1
        nums_array = evstr(rows(j))
        matrix = [matrix; nums_array]
    end
  end
endfunction

function res=Fnorm(m)
	res = 0.0
	for i = 1:size(m)(1)
		for j = 1:length(m(i))
			res = res + m(i,j) * m(i,j)
    end
  end
  res = sqrt(res)
endfunction

sizes = [3 100 200 300 400 500];
times = [];
absErrors = [];
relErrors = [];

[absErrorsQ, absErrorsR] = ([], []);
[relErrorsQ, relErrorsR] = ([], []);
for i = 1:length(sizes)
  if sizes(i) == 3 then
    A = readMatrixFromFile("matrixA_" + string(sizes(i)) + "_1.txt");
    myQ = readMatrixFromFile("matrixQ_" + string(sizes(i)) + "_1.txt");
    myR = readMatrixFromFile("matrixR_" + string(sizes(i)) + "_1.txt");
  else
    A = readMatrixFromFile("matrixA_" + string(sizes(i)) + ".txt");
    myQ = readMatrixFromFile("matrixQ_" + string(sizes(i)) + ".txt");
    myR = readMatrixFromFile("matrixR_" + string(sizes(i)) + ".txt");
  end
  tic()
  [Q, R] = qr(A);
  elapsed_time = toc();
  times = [times, elapsed_time];
  absErrors = [absErrors, Fnorm(Q * R - A)];
  relErrors = [relErrors, Fnorm(Q * R - A) / Fnorm(A)];
  
  absErrorsQ = [absErrorsQ, Fnorm(myQ - Q)];
  absErrorsR = [absErrorsR, Fnorm(myR - R)];
  relErrorsQ = [relErrorsQ, Fnorm(myQ - Q) / Fnorm(Q)];
  relErrorsR = [relErrorsR, Fnorm(myR - R) / Fnorm(R)];
end
for i = 1:length(sizes)
  disp("test size: " + string(sizes(i)))
  disp("time: " + string(times(i)))
  disp("abs error: " + string(absErrors(i)))
  disp("rel error: " + string(relErrors(i)))

  disp("abs error myQ: " + string(absErrorsQ(i)))
  disp("rel error myQ: " + string(relErrorsQ(i)))
 
  disp("abs error myR: " + string(absErrorsR(i)))
  disp("rel error myR: " + string(relErrorsR(i)))
  disp("-------------------------------------")
end
disp("all tests passed")
