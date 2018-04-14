function negD = HardestNeg(filepath, W, b, threshold)
    image_struct = dir(strcat(filepath, '/*.jpg'));      
    num_images = length(image_struct);
    load(strcat('/home/shivang/Desktop/CSE_512/hw2_q4/hw2data', '/trainAnno.mat'), 'ubAnno'); 
    
    ScoreValues = [];
    for i = 1:num_images
        imgs{i} = imread(strcat(strcat(filepath, '/'), image_struct(i).name));
        rects = HW2_Utils.detect(imgs{i}, W,b,0);
        rects(1:4,:) = cast(rects(1:4,:), 'uint8');
        I = rects(5,:) > 0;    
        rects(:, I==0) = [];
        
        ubs = ubAnno{i};
        for j=1:size(ubs,2)
            overlap = HW2_Utils.rectOverlap(rects, ubs(:,j));                    
            rects = rects(:, overlap < threshold);
            if isempty(rects)
                break;
            end
        end
        ScoreValues = cat(2, ScoreValues, rects(5, :));
%         if size(rects, 2) > 10
%             %flipped_rects =  fliplr(rects);
%             %rects = flipped_rects(:, 1:10);
%             rects = rects(:, 1:10);
%         end
        nNeg2SamplePerIm = size(rects, 2);
        im = imgs{i};
        [D_i, R_i] = deal(cell(1, nNeg2SamplePerIm));
        for j=1:nNeg2SamplePerIm
            imReg = im(rects(2,j):rects(4,j), rects(1,j):rects(3,j),:);
            imReg = imresize(imReg, HW2_Utils.normImSz);
            R_i{j} = imReg;
            D_i{j} = HW2_Utils.cmpFeat(rgb2gray(imReg));                    
        end
        
        negD{i} = cat(2, D_i{:});                
        negRegs{i} = cat(4, R_i{:});        
    end
    
    negD = cat(2, negD{:});
    if size(ScoreValues, 2) > 1000
        [m, i] = maxk(ScoreValues, 1000);
        negD = negD(:,i);
    end
    
    
%     size(negD)  
end

% % % % function [D, lb, imRegs] = getPosAndRandomNegHelper(dataset)
% % % %             rng(1234); % reset random generator. Keep same seed for repeatability
% % % %             load(sprintf('%s/%sAnno.mat', HW2_Utils.dataDir, dataset), 'ubAnno');
% % % %             [posD, negD, posRegs, negRegs] = deal(cell(1, length(ubAnno)));            
% % % %             
% % % %             for i=1:length(ubAnno)
% % % %                 ml_progressBar(i, length(ubAnno), 'Processing image');
% % % %                 im = imread(sprintf('%s/%sIms/%04d.jpg', HW2_Utils.dataDir, dataset, i));
% % % %                 %im = rgb2gray(im);
% % % %                 ubs = ubAnno{i}; % annotated upper body
% % % %                 if ~isempty(ubs)
% % % %                     [D_i, R_i] = deal(cell(1, size(ubs,2)));
% % % %                     for j=1:length(D_i)
% % % %                         ub = ubs(:,j);
% % % %                         imReg = im(ub(2):ub(4), ub(1):ub(3),:);
% % % %                         imReg = imresize(imReg, HW2_Utils.normImSz);
% % % %                         D_i{j} = HW2_Utils.cmpFeat(rgb2gray(imReg));
% % % %                         R_i{j} = imReg;
% % % %                     end 
% % % %                     posD{i}    = cat(2, D_i{:});                    
% % % %                     posRegs{i} = cat(4, R_i{:});
% % % %                 end
% % % %                 
% % % %                 % sample k random patches; some will be used as negative exampels
% % % %                 % Choose k sufficiently large to ensure success
% % % %                 k = 1000;
% % % %                 [imH, imW,~] = size(im);
% % % %                 randLeft = randi(imW, [1, k]);
% % % %                 randTop = randi(imH, [1, k]);
% % % %                 randSz = randi(min(imH, imW), [1, k]);
% % % %                 randRects = [randLeft; randTop; randLeft + randSz - 1; randTop + randSz - 1];
% % % %                 
% % % %                 % remove random rects that do not lie within image boundaries
% % % %                 badIdxs = or(randRects(3,:) > imW, randRects(4,:) > imH);
% % % %                 randRects = randRects(:,~badIdxs);
% % % %                 
% % % %                 % Remove random rects that overlap more than 30% with an annotated upper body
% % % %                 for j=1:size(ubs,2)
% % % %                     overlap = HW2_Utils.rectOverlap(randRects, ubs(:,j));                    
% % % %                     randRects = randRects(:, overlap < 0.3);
% % % %                     if isempty(randRects)
% % % %                         break;
% % % %                     end;
% % % %                 end;
% % % %                 
% % % %                 % Now extract features for some few random patches
% % % %                 nNeg2SamplePerIm = 2;
% % % %                 [D_i, R_i] = deal(cell(1, nNeg2SamplePerIm));
% % % %                 for j=1:nNeg2SamplePerIm
% % % %                     imReg = im(randRects(2,j):randRects(4,j), randRects(1,j):randRects(3,j),:);
% % % %                     imReg = imresize(imReg, HW2_Utils.normImSz);
% % % %                     R_i{j} = imReg;
% % % %                     D_i{j} = HW2_Utils.cmpFeat(rgb2gray(imReg));                    
% % % %                 end
% % % %                 negD{i} = cat(2, D_i{:});                
% % % %                 negRegs{i} = cat(4, R_i{:});
% % % %             end    
% % % %             posD = cat(2, posD{:});
% % % %             negD = cat(2, negD{:});   
% % % %             D = cat(2, posD, negD);
% % % %             lb = [ones(size(posD,2),1); -ones(size(negD,2), 1)];
% % % %             imRegs = cat(4, posRegs{:}, negRegs{:});            
% % % % end