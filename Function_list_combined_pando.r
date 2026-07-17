

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Base 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


### remove assay data 
###
### @input : obj , useAssay ="RNA",feature
### @return : None

resetAssayt <- function(obj , useAssay ="RNA",feature =NULL ){
  DefaultAssay(obj) <- useAssay
  assay_list <- obj@assays %>% names()
  print(assay_list)
  for ( a in assay_list){
    if ( a != as.character(useAssay) ){
      obj[[a]] <- NULL
    }
  }
  if( !is.null(feature) ) {
    obj = subset(obj, features = feature)
  }
  return(obj)
}

### Get obj list for clean obj
###
### @input : standardt_obj_list ,  obj
### @return : vars


GetRmObjList<- function( objs =NULL ,standardt_obj_list = NULL ){
  if( is.null(standardt_obj_list ) ){
    standardt_obj_list <- c("standardt_obj_list","param","SeuratOBJ","cluster_pallete", 
                            "palettedt","Category_list","subsetOBJ","subset_Cluster", "LINGER_data_list")
  }
  # すべてのオブジェクトを取得
  # 関数でないものだけ抽出
  vars <- objs[!sapply(objs, function(x) is.function(get(x)))]
  vars <- vars[ !vars %in% standardt_obj_list]
  # 結果を表示
  return(vars)
}



### ggplot png & pdf
###
### @input : p, outfile  height = 10,width = 10
### @return : None
wggplot <- function(p , file , height = 10,width = 10 ){
  outfile = file
  output.file <- paste(outfile ,".pdf", sep = "")
  ggsave( p , file =output.file , height = height,width = width )
  output.file <- paste(outfile ,".png", sep = "")
  ggsave( p , file =output.file , height = height,width = width )
}

### ggplotmake_outdir png & pdf
###
### @input : p, outfile  height = 10,width = 10
### @return : None

make_outdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("ディレクトリを作成しました: ", path)
  } else {
    message("output ディレクトリは既に存在します: ", path)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# EXP 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### Get Expression Rate matrix
##

GetExpRate <- function(obj , assay="RNA" ,layer="counts", cellLabel  ){
  use_cell <- intersect(names(cellLabel)  , colnames(obj) )%>% unlist
  cellLabel <- cellLabel[use_cell]
  expmat <- GetAssayData(obj[,use_cell] , assay = assay, layer = layer ) %>% 
    as.data.frame
  #expmat<- expmat[,expmat]
  categly_list <- unique(cellLabel)
  
  apply( expmat, 1, function(x , use_cell , categly_list){
    sapply(categly_list, function(y,x){
      extract_cell = names(cellLabel[cellLabel == y])
      pct = 100 * sum(x[extract_cell] > 0) /  length(extract_cell)
      return(pct)
    },x=x) %>% unlist
  } ,use_cell = use_cell, categly_list = categly_list ) %>%
    return()
}

calc_exprateByGroup <- function(rna_counts, cellLabel) {
  # rna_counts: genes x cells (dgCMatrix 推奨)
  if (!inherits(rna_counts, "dgCMatrix")) rna_counts <- as(rna_counts, "dgCMatrix")
  
  use_cell <- intersect(names(cellLabel)  , colnames(rna_counts) )%>% unlist
  cellLabel <- cellLabel[use_cell]
  rna_counts <-rna_counts[, use_cell ]
  categly_list <- unique(cellLabel)
  
  # >0 判定の疎行列を0/1化して rowSums
  bin <- rna_counts
  bin@x <- rep.int(1, length(bin@x))
  
  sapply( categly_list , function(x, rna_counts ,use_cell){
    extract_cell = names(use_cell[use_cell == x])
    pct = Matrix::rowSums(bin[, extract_cell, drop = FALSE]) /  length(extract_cell)
  },categly_list = categly_list ,use_cell=use_cell )%>%
    return()
}  
  
#  apply( bin , 1, function(x , use_cell , categly_list){
#    sapply(categly_list, function(y,x){
#      extract_cell = names(cellLabel[cellLabel == y])
#      pct = Matrix::rowSums(bin[, extract_cell, drop = FALSE]) /  length(extract_cell)
#      return(pct)
#    },x=x) %>% unlist
#  } ,use_cell = use_cell, categly_list = categly_list ) %>%
#    return()

  


# RNA counts と cell名から、human/marmoset の発現率をベクトル化で計算
calc_exprate <- function(rna_counts, human_pattern = "Human") {
  # rna_counts: genes x cells (dgCMatrix 推奨)
  if (!inherits(rna_counts, "dgCMatrix")) rna_counts <- as(rna_counts, "dgCMatrix")
  cell_is_human <- grepl(human_pattern, colnames(rna_counts), fixed = TRUE)
  nh <- sum(cell_is_human)
  nm <- sum(!cell_is_human)
  
  # >0 判定の疎行列を0/1化して rowSums
  bin <- rna_counts
  bin@x <- rep.int(1, length(bin@x))
  
  human_open <- Matrix::rowSums(bin[, cell_is_human, drop = FALSE])
  marm_open  <- Matrix::rowSums(bin[, !cell_is_human, drop = FALSE])
  
  data.frame(
    human = human_open / nh,
    marmoset = marm_open / nm,
    row.names = rownames(rna_counts)
  )
}

### MakeCorrelationPlot
###
### @input :     Compare_obj_list ,outdir ,common_gene , category 
### @return : None

MakeCorrelationPlot <- function( Compare_obj_list ,outdir ,common_gene , category = "bin" ,assayD="RNA"){
  if (  !category %in% colnames(Compare_obj_list[[1]]@meta.data) ) {
    print( paste0( "category: ", category , " not inmetadate. plese retry..") )
    return(0)
  }
  ## main comapre 2 species
  # resset assay
  Compare_obj_list <- lapply(Compare_obj_list , function(x){ resetAssayt(x, assayD ,common_gene) } )
  spe_name = names(Compare_obj_list)
  
  ## 平均発現量での相関を評価
  mean_bin_exp_spe1 = AverageExpression( Compare_obj_list[[1]] , assays = assayD ,group.by =  category )[[assayD]]
  mean_bin_exp_spe2 = AverageExpression( Compare_obj_list[[2]] , assays = assayD ,group.by =  category )[[assayD]]
  
  colnames(mean_bin_exp_spe1) <- paste0( gsub("[0-9]$","",spe_name[1]) , "_" , gsub("^g","",colnames(mean_bin_exp_spe1)) )
  colnames(mean_bin_exp_spe2) <- paste0( gsub("[0-9]$","",spe_name[2])  , "_" , gsub("^g","",colnames(mean_bin_exp_spe2)) )
  bin_cor_matrix <- cor(mean_bin_exp_spe2[common_gene,] %>% as.matrix() , mean_bin_exp_spe1[common_gene,] %>% as.matrix() )
  
  # barplot of bin count 
  if ( category == "bin" | category == "mod_bin") {
    
    if (category == "bin" ){
      bin_count_spe1 <-table(Compare_obj_list[[1]]@meta.data$RNA_cluster, 
                             Compare_obj_list[[1]]@meta.data$bin) %>% as.matrix() %>% t
      
      bin_count_spe2 <-table(Compare_obj_list[[2]]@meta.data$RNA_cluster, 
                             Compare_obj_list[[2]]@meta.data$bin) %>% as.matrix() %>% t
      rownames(bin_count_spe1) <- 0:43    
      rownames(bin_count_spe2) <- 0:43
    } else{
      
      bin_count_spe1 <-table(Compare_obj_list[[1]]@meta.data$RNA_cluster, 
                             Compare_obj_list[[1]]@meta.data$mod_bin) %>% as.matrix() %>% t
      
      bin_count_spe2 <-table(Compare_obj_list[[2]]@meta.data$RNA_cluster, 
                             Compare_obj_list[[2]]@meta.data$mod_bin) %>% as.matrix() %>% t
      print( bin_count_spe1)
      rownames(bin_count_spe1)[1:44] <- 0:43
      rownames(bin_count_spe2)[1:44]<- 0:43
    }
    #+++ Heatamap ++++++
    ## col annotaion = x軸 spe2
    #HeatmapAnnotaion 
    annotation_spe1 <-  setNames( 
      list( anno_barplot(bin_count_spe1,bar_width = 1,
                         gp = gpar(fill = cluster_pallete[colnames(bin_count_spe1)] ),
                         border = TRUE,axis = TRUE,
                         height = unit(4, "cm"))  ) 
      , spe_name[1] )  
    annotation_spe2 <- setNames( 
      list( anno_barplot(bin_count_spe2,bar_width = 1,
                         gp = gpar(fill = cluster_pallete[colnames(bin_count_spe2)] ),
                         border = TRUE,axis = TRUE,which = "row",
                         width = unit(4, "cm"))  ) 
      , spe_name[2] )  
  } else {
    bin_count_spe1 <- table(Compare_obj_list[[1]]@meta.data[, category])  %>% as.matrix
    #rownames(bin_count_spe1) <- 0:43
    bin_count_spe2 <-table(Compare_obj_list[[2]]@meta.data[, category])  %>% as.matrix
    #rownames(bin_count_spe2) <- 0:43
    
    annotation_spe1 <-  setNames( 
      list( anno_barplot(bin_count_spe1,bar_width = 1,
                         gp = gpar(fill = cluster_pallete[rownames(bin_count_spe1)] ),
                         border = TRUE,axis = TRUE,
                         height = unit(4, "cm"))  ) 
      , spe_name[1] )  
    annotation_spe2 <- setNames( 
      list( anno_barplot(bin_count_spe2,bar_width = 1,
                         gp = gpar(fill = cluster_pallete[rownames(bin_count_spe2)] ),
                         border = TRUE,axis = TRUE,which = "row",
                         width = unit(4, "cm"))  ) 
      , spe_name[2] )  
  }
  
  ha <- do.call(HeatmapAnnotation, annotation_spe1)
  row_ha <- do.call(rowAnnotation ,annotation_spe2 )
  
  col_viridis = viridis(n = 44, option = "turbo") 
  
  #print(bin_cor_matrix )
  ht=Heatmap(bin_cor_matrix,
             name = "Correlation", 
             col = col_viridis ,
             top_annotation = ha, 
             left_annotation = row_ha,
             show_column_names = T,
             show_row_names =  T,
             cluster_columns = FALSE,
             cluster_rows = FALSE , width = unit(18, "cm"),height = unit(18, "cm"))
  #+++ Heatamap ++++++
  
  output.file = paste0( outdir,"bincorrelationHeatmap.pdf" )
  pdf( file = output.file , height = 14 ,width = 14)
  
  draw(ht)
  dev.off()
  
  output.file = paste0( outdir,"bincorrelationHeatmap.png" )
  png( file = output.file , height = 1000 ,width = 1000)
  draw(ht)
  dev.off()
  return( list( cormat = bin_cor_matrix , mean1 = mean_bin_exp_spe1 , mean2 = mean_bin_exp_spe2  ) )
}


### 転写開始点周辺の領域を定義する関数 , 上流　下流を自由に定義可能
###
### @input : param,SeuratOBJ , Category_list
### @return : SeuratOBJ 
resize_by_strand_with_region_string <- function(gr, up = 0, down = 0) {
  #half_width <- width / 2
  
  # 元の start/end を保存
  original_start <- start(gr)
  original_end <- end(gr)
  new_start <- ifelse(strand(gr) == "+",
                      original_start - up,
                      original_end - down -1)
  new_end <- ifelse(strand(gr) == "+",
                    original_start + down +1,
                    original_end + up )
  # 新しい GRanges を作成
  gr_new <- GRanges(seqnames = seqnames(gr),
                    ranges = IRanges(start = new_start, end = new_end),
                    strand = strand(gr))
  
  # 元のメタデータをコピー
  mcols(gr_new) <- mcols(gr)
  
  # 元の領域情報を "chr1-718463-719301" のような文字列で追加
  region_string <- paste0(as.character(seqnames(gr)), "-", original_start, "-", original_end)
  mcols(gr_new)$original_region <- region_string
  
  return(gr_new)
}



###　Convert Fitted gene data to matrix
###
### @input : param,SeuratOBJ , Category_list
### @return : SeuratOBJ 
converSplineToMat <- function( spline.fits , bin_range= 0:43  ,gene_col = "GeneID"){
  #spline.fits <- spline.fits %>% mutate(x = as.numeric(x) , y=as.numeric(y))
  spline.fits.mat <-  sapply(spline.fits ,function(x){
    #x <- x %>% mutate( x = as.numeric(x) , y=as.numeric(y)) 
    pivot_wider(x[,2:3] ,names_from = 'x', values_from = 'y' ) %>%
      as.numeric() })  %>% as.data.frame() %>% t() %>% as.data.frame()
  
  if ( length(bin_range) != ncol(spline.fits.mat) ) {
    print(" bin_rangeg is diffrent from matrix size.")
    return(NULL)
  }
  
  rownames(spline.fits.mat) <-sapply(spline.fits,function(x){  x[,gene_col][1] })
  colnames(spline.fits.mat) <- bin_range
  spline.fits.mat.scale <- spline.fits.mat %>% t %>% scale()%>%t%>% as.data.frame()
  #print(  spline.fits.mat.scale %>% head )
  maxbin <- apply(spline.fits.mat.scale,1,function(x){ which.max(x)})%>% as.numeric() -1 
  
  return( list( mat = spline.fits.mat, scalemat = spline.fits.mat.scale, maxbin = maxbin)  )
}





#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# DEG
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

make_gene_by_cluster_status <- function(df,
                                        use_padj = TRUE,
                                        padj_threshold = 0.05) {
  # cluster を文字化（ワイド列名に使うため）
  df2 <- df %>%mutate(cluster_chr = as.character(cluster))
  # Cluster一覧（欠番なく並べたい場合は sort）
  clusters <- df2 %>%
    distinct(cluster_chr) %>%
    pull(cluster_chr) %>%
    sort()
  
  # （任意）有意フィルタ：検出の定義を p_val_adj で絞る
  if (use_padj) {
    df2 <- df2 %>%
      filter(!is.na(p_val_adj), p_val_adj < padj_threshold)
  }
  
  # gene × cluster のステータスを作成
  status_long <- df2 %>%
    mutate(status = case_when(
      avg_log2FC > 0 ~ "UP",
      avg_log2FC < 0 ~ "DOWN",
      TRUE ~ NA_character_   # 0やNAは NA とする
    )) %>%
    dplyr::select(gene, cluster_chr, status) %>%
    # 未検出（行がない）組み合わせは NA で補完
    complete(gene, cluster_chr = clusters, fill = list(status = NA_character_))
  
  # ワイド化：Cluster0, Cluster1, ... 列にステータス
  status_wide <- status_long %>%
    pivot_wider(
      names_from  = cluster_chr,
      values_from = status,
      names_prefix = "Cluster"
    ) %>%
    arrange(gene)
  
  joned_status_wide <- df %>%
    left_join(status_wide, by = "gene")
  
  return(joned_status_wide)
}

add_otheceClusterDEG <- function( deglist ,use_col ="cluster" ){
  deglist <- deglist %>%
    group_by(gene) %>%
    mutate(UpCluster =.data[[use_col]][avg_log2FC > 0] %>% unique()  %>% sort %>% paste(.,collapse = ","),
           DownCluster =.data[[use_col]][avg_log2FC < 0] %>% unique()  %>% sort %>% paste(.,collapse = ",") ) %>%
    ungroup()
  return(deglist)
}




#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Pseudtime
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

###  GAM fitting 
###
### @input : param,SeuratOBJ , Category_list
### @return : SeuratOBJ 

fitSmoothGAM  <- function(Y , genes , cell.ord , pseudotime ,num.cores =5 ,L1,total_day = 42 ,bin_day= 1,k =6){
  library(mgcv)
  expr_matrix <-Y[genes, cell.ord]
  sc.tc.df.adj <- expr_matrix  %>% 
    t %>% 
    as.data.frame() %>% 
    rownames_to_column( var ="Cell") %>% 
    pivot_longer(-Cell,names_to = "GeneID",values_to = "log2.expr") %>% 
    left_join( . , (L1$sds.data %>% dplyr::select( Cell , t ,rep)) , by = "Cell" ) %>%
    #mutate( t= (t -min(t))/ ( max(t)-min(t))) %>%
    transmute(y = log2.expr, tme = t, ind = rep, variable = GeneID)  
  
  spline.fits <- mclapply(1:length(genes), function(i){
    # GAMで周期的スプラインをフィッティング
    tmp <-  sc.tc.df.adj %>% dplyr::filter(variable == genes[i]) %>%
      transmute(GeneID = variable, x = tme, y = y)
    y <- tmp$y
    t <- tmp$x
    fit <- gam(y ~ s(t, bs = "cc",k = k))  
    newdata <- data.frame(t = seq(0, total_day, by = bin_day))
    rna.sp <- predict(fit, newdata = newdata)
    mu <- data.frame(x = seq(0, total_day, by = bin_day), y = rna.sp) 
    mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
    return(mu)
  }, mc.cores = num.cores)
  
  return(spline.fits)
}

## Input  : param ,subsetOBJ  
## Output : timeBin_data_org

### pseudotime Calculation 
###
### @input : subsetOBJ,reverset.time = T ,total_day = 42 ,bin_day= 1,start_cells
### @return : list(  sds.data ,umap.sds )
compute_pseudotime_on_ellipse <- function(
    seu,
    ell,
    umap_reduction = "umap",
    start_cell = NULL,
    use_arclength = FALSE,
    clockwise = FALSE
) {
  # --- UMAP 埋め込みの取り出し ---
  umap_df <- Embeddings(seu, reduction = umap_reduction) %>%
    as.data.frame() %>%
    rownames_to_column("Cell")
  # UMAP 列名（UMAP_1/UMAP_2 想定。異なる場合はここを調整）
  umap_cols <- intersect(c("UMAP_1", "UMAP_2", "umap_1", "umap_2", "U1", "U2", "x", "y"),
                         colnames(umap_df))
  stopifnot(length(umap_cols) >= 2)
  umap_cols <- umap_cols[1:2]
  
  # --- 4) 楕円の中心座標 ---
  # 1) Ellipsefit の出力の中心を使用
  center <-  c(  ell$Para$Centerpoint_X, ell$Para$Centerpoint_Y)
  names(center) <- umap_cols
  
  # --- 5) 各細胞の角度（パラメトリック角 t） ---
  #   A) 簡易：回転や長短軸を無視して「中心から見た極角」
  #   B) 精密：楕円の向きと長短軸（a,b）を推定し、t = atan2( (y')/b, (x')/a )
  #      （ell に角度や a,b があればそれを使い、無ければ PCA で推定）
  
  # a,b,phi を PCA で推定
  phi <- NULL; a <- NULL; b <- NULL
  coord <- as.data.frame(ell$Coord)
  coord_cols <- intersect(umap_cols, colnames(coord))
  if (length(coord_cols) < 2) coord_cols <- colnames(coord)[1:2]
  M <- sweep(as.matrix(coord[, coord_cols, drop = FALSE]), 2, center, "-")
  pc <- prcomp(M, center = FALSE, scale. = FALSE)
  # 第一主成分の向き
  phi <- atan2(pc$rotation[2, 1], pc$rotation[1, 1])
  # 回転を元に、軸方向での最大半径を a,b とする（外接長半径推定）
  Rm <- matrix(c(cos(-phi), -sin(-phi), sin(-phi), cos(-phi)), 2, 2, byrow = TRUE)
  Mr <- M %*% t(Rm)
  a <- max(abs(Mr[, 1]))  # 長軸半径
  b <- max(abs(Mr[, 2]))  # 短軸半径
  
  # 細胞ごとのベクトル
  V <- as.matrix(umap_df[, umap_cols]) - matrix(center, nrow(umap_df), 2, byrow = TRUE)
  
  # パラメトリック角 t を計算
  # 角度のみで良い場合：phi, a, b を使って回転・縮尺を補正した atan2((y')/b, (x')/a)
  Rm <- matrix(c(cos(-phi), -sin(-phi), sin(-phi), cos(-phi)), 2, 2, byrow = TRUE)
  Vr <- V %*% t(Rm)
  t_param <- atan2(Vr[, 2] / b, Vr[, 1] / a)
  t_param <- (t_param + 2 * pi) %% (2 * pi)  # 0〜2π
  
  # 擬時間（角度ベース）
  pt_angle <- t_param / (2 * pi)  # 0〜1
  
  # （任意）弧長正規化：0〜t の弧長 / 全弧長
  if (use_arclength) {
    integrand <- function(u) sqrt((a * sin(u))^2 + (b * cos(u))^2)
    # 総弧長を1回だけ積分
    total_len <- integrate(integrand, lower = 0, upper = 2 * pi,
                           subdivisions = 400L, rel.tol = 1e-6)$value
    arc_len <- map_dbl(t_param, ~ integrate(integrand, lower = 0, upper = .x,
                                            subdivisions = 200L, rel.tol = 1e-6)$value)
    pt <- arc_len / total_len
  } else {
    pt <- pt_angle
  }
  
  # --- 6) スタート細胞で 0 にリセット（向き反転も可） ---
  # start_cell が無ければ「中心から見て +x 方向（t ≈ 0）に最も近い細胞」を採用
  if (is.null(start_cell)) {
    start_idx <- which.min(abs(t_param - 0))
  } else {
    stopifnot(start_cell %in% umap_df$Cell)
    start_idx <- which(umap_df$Cell == start_cell)[1]
  }
  t0 <- pt[start_idx]
  pt_shifted <- (pt - t0) %% 1
  if (clockwise) {
    pt_shifted <- (1 - pt_shifted) %% 1
  }
  
  # 結果を Seurat に格納
  seu@meta.data$pseudotime_ellipse <- pt_shifted
  seu@meta.data$pseudotime_ellipse_angle <- pt_angle # 角度ベースの素の値も残す
  
  list(
    seurat = seu,
    center = setNames(center, umap_cols),
    a = a, b = b, phi = phi,
    t_param = t_param,
    pseudotime = pt_shifted,
    df = dplyr::bind_cols(umap_df, tibble::tibble(
      t_param = t_param,
      pseudotime_ellipse = pt_shifted
    )),
    start_cell = umap_df$Cell[start_idx]
  )
}

#!
#可視化（擬時間で着色）
#!
plot_pseudotime_umap <- function(df, umap_cols = c("UMAP_1","UMAP_2"),
                                 pt_col = "pseudotime_ellipse",
                                 center = NULL, start_cell = NULL) {
  p <- ggplot(df, aes_string(umap_cols[1], umap_cols[2], color = pt_col)) +
    geom_point(size = 0.8, alpha = 0.8) +
    scale_color_viridis_c(option = "turbo") +
    coord_equal() +
    theme_classic()
  if (!is.null(center)) {
    p <- p + 
      geom_point(aes(x = center[1], y = center[2]),
                 inherit.aes = FALSE, shape = 4, size = 3, color = "black", stroke = 1.2)
  }
  if (!is.null(start_cell) && start_cell %in% df$Cell) {
    sc <- df[df$Cell == start_cell, umap_cols]
    p <- p + 
      geom_point(aes(x = sc[[1]], y = sc[[2]]),
                 inherit.aes = FALSE, shape = 21, size = 2.4, stroke = 0.8,
                 fill = "white", color = "black")
  }
  p
}



fitPseudoTime_theta <- function(subsetOBJ, 
                                reverset.time = T ,
                                total_day = 42 ,bin_day= 1,
                                start_cells = "I7704M_1yr8mnth_CCTAAGTAGGCTGGCT-1"){
  #++++++++++++++++++++++++
  # Fitting curve
  #++++++++++++++++++++++++
  ## Sample == Cell 
  umap_emb <- subsetOBJ@reductions$umap@cell.embeddings %>% as.data.frame() %>% dplyr::mutate(Cell = rownames(.)) 
  umap_emb$cluster <- subsetOBJ$seurat_clusters
  ell <- Ellipsefit(umap_emb,umap_1, umap_2, coords = TRUE, bbox = TRUE)
  coord <- ell$Coord
  
  
  x <- as.matrix(umap_emb[,c(1,2)])
  fit <- principal_curve(x, start = as.matrix(coord), smoother = "periodic_lowess", maxit = 0)
  tmp_pt <- fit$lambda
  
  
  # 角度ベース（軽量）で擬時間
  compute_pseudotime_res <- compute_pseudotime_on_ellipse(
    seu = subsetOBJ,
    ell = ell,
    umap_reduction = "umap",
    start_cell = start_cells,
    use_arclength = T,  # TRUE にすると弧長正規化
    clockwise = reverset.time
  )
  
  # 結果の取り出し
  subsetOBJ_T <- compute_pseudotime_res$seurat
  head(subsetOBJ_T@meta.data$pseudotime_ellipse)
  
  # 可視化
  p <- plot_pseudotime_umap(
    df = compute_pseudotime_res$df,
    umap_cols = names(compute_pseudotime_res$center),
    pt_col = "pseudotime_ellipse",
    center = compute_pseudotime_res$center,
    start_cell = compute_pseudotime_res$start_cell
  )
  
  output.file = paste0( param$outdir,"/Psedotime_tan_predection" )
  wggplot( p, file =output.file , height = 5,width = 5 )
  
  
  #++++++++++++++++++++++++
  ### change start cell
  #++++++++++++++++++++++++
  #start_pt <- tmp_pt[start_cells]
  #tmp_pt <- tmp_pt - start_pt
  #tmp_pt[ tmp_pt <0 ] <- tmp_pt[ tmp_pt <0 ] + max(abs(tmp_pt[ tmp_pt <0 ]))  + max(tmp_pt)
  #tmp_pt <- (tmp_pt - min(tmp_pt))/(max(tmp_pt) - min(tmp_pt))
  
  tmp_pt <-compute_pseudotime_res$df%>% column_to_rownames(var = "Cell") %>% dplyr::select(pseudotime_ellipse) 
  
  p1 <- cbind(umap_emb,tmp_pt = tmp_pt[rownames(umap_emb),])  %>% 
    ggplot() + geom_point(aes( x=umap_1 , y=umap_2 , col =tmp_pt ), size = 0.5) +
    geom_path( data = coord , aes(x = x, y = y),colour = "red") +
    theme_classic()
  
  p2 <- cbind(umap_emb,tmp_pt = tmp_pt[rownames(umap_emb),])  %>% 
    mutate(  St = if_else( Cell == start_cells, "Start","other" )) %>%
    ggplot() + geom_point(aes( x=umap_1 , y=umap_2 , col =St ), size = 0.5) +
    geom_path( data = coord , aes(x = x, y = y),colour = "red") +
    theme_classic()
  print(p1 + p2)
  
  cell.ord <- order(tmp_pt$pseudotime_ellipse)
  sds.data <- data.frame(Cell = rownames(fit$s), 
                         cell.ord = cell.ord,
                         pt = tmp_pt$pseudotime_ellipse )
  #sc1 = fit$s[,1], 
  #sc2 = fit$s[,2]) %>% head
  
  #++++++++++++++++++++++++
  # change time to 42 day
  #++++++++++++++++++++++++
  pt <- sds.data$pt
  #if(reverset.time){
  #  pt <- max(pt) - pt
  #}
  
  #total_day = 42 ,bin_day = 1
  #Map the pseudo-time to 0-total day: bin_day
  t <- (total_day + bin_day) * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))
  sds.data$t.org <- t
  
  ## time-index cells in 1 day intervals and identify cells in each partition
  ## They will be considered as replicates
  time.breaks <- seq(bin_day, total_day + bin_day, by = bin_day) 
  time.idx <- rep(0, nrow(sds.data))
  
  ind <- which(sds.data$t.org <= time.breaks[1])
  time.idx[ind] <- 0
  
  for(i in 2:length(time.breaks)){
    ind <- which(sds.data$t.org > time.breaks[(i-1)] & sds.data$t.org <= time.breaks[i])
    time.idx[ind] <- i - 1
  }
  
  sds.data$time.idx <- time.idx
  ## Update the time to 1 day min increments
  sds.data$t <- sds.data$t.org
  sds.data <- sds.data %>%  
    group_by(time.idx) %>% mutate(rep = seq(1:n()))
  
  cbind(umap_emb, t = sds.data[,"time.idx"])%>% 
    ggplot() + geom_point(aes( x=umap_1 , y=umap_2 , col =time.idx ), size = 0.5) +
    geom_path( data = coord , aes(x = x, y = y),colour = "red") +
    theme_classic()
  
  L <- list( 
    sds.data = sds.data,
    umap.sds  = left_join(umap_emb, sds.data, by = "Cell"))
  return(L) 
}




###
### @input : param, subsetOBJ, total_day = 43 , num.cores = 5
### @return : fitting = sc.spline.fits, Gene_stats=Gene_stats

DoExpgittingcurve_species <- function( param, subsetOBJ, total_day = 43 , #外挿分を+1する
                                       num.cores = 8 , assay ="RNA"){
  timeBin_data <-readRDS( param$outfiles$out_timebin )
  DefaultAssay(subsetOBJ) <- assay
  # Human
  human_obj <- subset(subsetOBJ , species == "human")
  if(  dim(SeuratObject::LayerData(human_obj[[assay]], layer = "data"))[1] == 0 ){
    if(  assay == "RNA" ){
      human_obj <- NormalizeData(human_obj )
    } else if (assay == "SCT"){
      human_obj <- SCTransform(human_obj)
    }
  }
  cellname.human <- colnames(human_obj) 
  timeBin_data.human <- list( 
    sds.data = timeBin_data$sds.data %>% dplyr::filter( Cell %in% cellname.human),
    umap.sds = timeBin_data$umap.sds %>% dplyr::filter( Cell %in% cellname.human)
  )
  print("human fitting...")
  result_fit.human <- DoExpgittingcurve( param, subsetOBJ=human_obj, total_day = total_day , #外挿分を+1する
                                         num.cores = num.cores , timeBin_data =timeBin_data.human )
  # Marmoset
  marmoset_obj <- subset(subsetOBJ , species == "marmoset")
  if(  dim(SeuratObject::LayerData(marmoset_obj[[assay]], layer = "data"))[1] == 0 ){
    if(  assay == "RNA" ){
      marmoset_obj <- NormalizeData(marmoset_obj )
    } else if (assay == "SCT"){
      marmoset_obj <- SCTransform(marmoset_obj)
    }
  }
  cellname.marmoset <- colnames(marmoset_obj) 
  timeBin_data.marmoset <- list( 
    sds.data = timeBin_data$sds.data %>% dplyr::filter( Cell %in% cellname.marmoset),
    umap.sds = timeBin_data$umap.sds %>% dplyr::filter( Cell %in% cellname.marmoset)
  )
  print("marmoset fitting...")
  result_fit.marmoset <- DoExpgittingcurve( param,subsetOBJ = marmoset_obj, total_day = total_day , #外挿分を+1する
                                            num.cores = num.cores , timeBin_data =timeBin_data.marmoset )  
  #Human_res <- list(  )
  #Marmoset_res <- list()
  return( list(marmoset = result_fit.marmoset , human = result_fit.human))
}


### DoExpgittingcurve 
###
### @input : param, subsetOBJ, total_day = 43 , num.cores = 5
### @return : fitting = sc.spline.fits, Gene_stats=Gene_stats

DoExpgittingcurve<- function( param, subsetOBJ, total_day = 43 , #外挿分を+1する
                              num.cores = 8 ,  timeBin_data = NULL){
  if( is.null(timeBin_data) ){
    timeBin_data <-readRDS( param$outfiles$out_timebin )
  }
  ## select variable gene make cycle matrix
  Y <- subsetOBJ@assays$RNA$data 
  vat.genes_check <-apply(Y, 1, function(x){
    return( c( var =var(x) , mean=mean(x), exppct = sum(x>0)/length(x) ) ) }) %>%
    t
  summary(vat.genes_check )
  
  #+ check Gene stats
  Gene_stats <- vat.genes_check %>% 
    as.data.frame() %>% 
    rownames_to_column(var="Gene") %>%
    pivot_longer( -Gene , names_to = "Type", values_to = "val" )
  
  p1 <- Gene_stats %>%
    filter(Type == "var") %>%
    ggplot() +
    geom_histogram(aes(x = val),  fill = "red", binwidth = 0.01) +
    theme_classic()+
    theme(text = element_text(size = 20))
  
  p2 <- Gene_stats %>%
    filter(Type == "mean") %>%
    ggplot() +
    geom_histogram(aes(x = val), fill = "blue", binwidth = 0.2) +
    #scale_y_log10() +
    theme_classic()+
    theme(text = element_text(size = 20))
  
  p3 <- Gene_stats %>%
    filter(Type == "exppct") %>%
    ggplot() +
    geom_histogram(aes(x = val), fill = "green", binwidth = 0.1) +
    #scale_y_log10() +
    theme_classic()+
    theme(text = element_text(size = 20))
  p<-p1 | p2 |p3 
  
  output.file = paste0( param$outdir,"/GenExpStat.pdf" )
  wggplot( p, file =output.file , height = 4,width = 13 )
  
  #+ check Gene stats
  # ほとんど発現していない遺伝子を除く
  var.genes <- vat.genes_check %>% 
    as.data.frame %>%
    filter( var > 0.01 & mean >0.2  & exppct > 0.1 ) %>% 
    rownames()
  
  print(summary(vat.genes_check[var.genes,] ))
  
  var.genes.gamp <- GeneSelection(subsetOBJ , var.genes , mc.cores =num.cores ,ptd= timeBin_data$sds.data$t)
  var.genes.sig <-var.genes[var.genes.gamp < 0.01]
  
  cell.ord <- timeBin_data$sds.data$Cell
  pseudotime <- timeBin_data$sds.data %>% pull(time.idx) 
  names(pseudotime) <- timeBin_data$sds.data %>% pull( Cell)
  
  sc.spline.fits.spline <- fitSmoothSpline(Y ,genes = var.genes.sig , 
                                           cell.ord ,
                                           pseudotime , 
                                           total_day = total_day , #外挿分を+1する
                                           num.cores = num.cores ,
                                           L1=timeBin_data)
  sc.spline.fits.gam <-  fitSmoothGAM(Y ,genes = var.genes.sig , 
                                      cell.ord ,
                                      pseudotime , 
                                      total_day = total_day , #外挿分を+1する
                                      num.cores = num.cores ,
                                      L1=timeBin_data , k = 7)
  
  #  output.file <- paste0( param$outfiles$out_fitgene)
  #  saveRDS(sc.spline.fits , output.file)
  
  return( list( seuratobj = subsetOBJ,
                time.bin.data = timeBin_data ,
                fitting = sc.spline.fits.spline , 
                fitting.gam = sc.spline.fits.gam , 
                Gene_stats=Gene_stats) )
}



### gam 
###
### @input : subsetOBJ , var.genes , mc.cores = 5,ptd
### @return : gam.pval.adj 
GeneSelection <- function(subsetOBJ , var.genes , mc.cores = 5,ptd){
  #　新しく遺伝子をセレクション
  ## select variable gene make cycle matrix
  Y <- subsetOBJ@assays$RNA$data 
  Y <- Y[var.genes, ]
  print(Y %>%dim)
  cat('Fitting the GAM model\n')
  gam.pval <- mclapply(1:nrow(Y), function(z){
    d <- data.frame(z=as.numeric(Y[z,]), t=as.numeric(ptd))
    tmp <- gam::gam(z ~ lo(t), data=d)
    p <- summary(tmp)$anova$`Pr(F)`[2] ## nonlinear effects
    p
  }, mc.cores = mc.cores)
  
  gam.pval <- unlist(gam.pval)
  #print(gam.pval )
  names(gam.pval) <- rownames(Y)
  ## Remove the NA's and get the best fits
  if(any(is.na(gam.pval))){
    gam.pval <- gam.pval[-which(is.na(gam.pval))]
  }
  
  gam.pval.adj <- p.adjust(gam.pval, method = 'fdr', n = length(gam.pval))
  gam.pval.sig <- gam.pval[gam.pval.adj < 0.01] 
  print(length(gam.pval.sig)) ## number of correlating genes
  
  return(gam.pval.adj)
}



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# PANDO 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

### pand function
### Pando GRN - カスタムノードカラー UMAP 可視化
### 依存: ggraph, igraph, tidygraph, ggplot2, dplyr


###' Pando GRN を任意のノードカラーで UMAP 可視化
###'
###' @param grn         Pando の GRN オブジェクト (Seurat ベース)
###' @param netpath     NetworkGraph に渡すグラフ名 (デフォルト: "umap_graph")
###' @param node_values named numeric vector  例: c(GATA1=1.5, TAL1=-0.3)
###'                    NULL の場合は degree centrality を使用（デフォルト）
###' @param node_color  カテゴリカラー用 named character vector
###'                    例: c(GATA1="TF", TAL1="Target")
###'                    node_values と同時指定した場合は node_values が優先
###' @param palette     連続値: RColorBrewer / viridis パレット名 (例: "viridis", "RdBu")
###'                    カテゴリ: 色ベクトル 例: c(TF="#E64B35", Target="#4DBBD5")
###' @param color_label カラーバー / 凡例のラベル
###' @param node_size         ノードサイズ (centrality によるサイズ変更をしない場合の固定値)
###' @param size_range        centrality によるサイズ範囲 c(最小, 最大)  例: c(1, 8)
###' @param edge_alpha        エッジの透明度
###' @param edge_width_range  estimate絶対値によるエッジ幅の範囲 c(最小, 最大)  例: c(0.2, 2)
###' @param title       図タイトル
###' @param show_labels "none"  : ラベルなし（デフォルト）
###'                    "all"   : 全ノードにラベル表示
###'                    character vector: 表示したい遺伝子名を指定
###'                    例: c("GATA1", "TAL1", "MYC")
###' @param label_size  ラベルのフォントサイズ (デフォルト: 3)
###' @param label_color ラベルの文字色 (デフォルト: "black")
###'
###' @return ggplot オブジェクト

plot_grn_custom <- function(
    grn,
    netpath      = "umap_graph",
    node_values  = NULL,   # named numeric: 連続値カラー
    node_color   = NULL,   # named character: カテゴリカラー
    palette      = "viridis",
    color_label  = "Value",
    node_size    = 3,
    size_range        = c(1, 8),    # centrality によるノードサイズ範囲
    edge_alpha        = 0.6,
    edge_width_range  = c(0.2, 2),  # estimate絶対値によるエッジ幅の範囲
    title        = "GRN UMAP",
    show_labels  = "none",  # "none" / "all" / character vector of gene names
    label_size   = 3,
    label_color  = "black"
) {
  
  # ── 1. Pando から igraph を取り出す ──────────────────────────────────────
  g <- NetworkGraph(grn, graph = netpath)
  
  # ── 2. UMAP 座標の取得 (vertex attribute から) ───────────────────────────
  umap_coords <- data.frame(
    UMAP1 = as.vector(V(g)$UMAP_1),
    UMAP2 = as.vector(V(g)$UMAP_2),
    row.names = V(g)$name
  )
  lay <- as.matrix(umap_coords[, 1:2])
  
  # ── 3. ノードカラー値・サイズ・ラベルを igraph に付与 ───────────────────
  genes <- V(g)$name
  
  # centrality でノードサイズをスケーリング
  cent <- V(g)$centrality
  cent_scaled <- scales::rescale(cent, to = size_range)
  V(g)$node_size <- cent_scaled
  
  if (!is.null(node_values)) {
    V(g)$color_val <- node_values[genes]
    color_type <- "continuous"
  } else if (!is.null(node_color)) {
    V(g)$color_cat <- node_color[genes]
    color_type <- "categorical"
  } else {
    # デフォルト: degree centrality
    V(g)$color_val <- degree(g, normalized = TRUE)
    color_type <- "continuous"
    color_label <- "Degree Centrality"
  }
  
  # ラベル: 表示対象のみ遺伝子名、それ以外は NA
  if (identical(show_labels, "none")) {
    V(g)$label_name <- NA_character_
  } else if (identical(show_labels, "all")) {
    V(g)$label_name <- genes
  } else {
    V(g)$label_name <- ifelse(genes %in% show_labels, genes, NA_character_)
  }
  
  # ── 4. エッジの estimate をグラフに付与 ─────────────────────────────────
  # estimate の絶対値を 0.2〜2 の幅にスケーリング
  est <- E(g)$estimate
  est_width <- scales::rescale(abs(est), to = edge_width_range)
  E(g)$est_width <- est_width
  E(g)$est_sign <- est   # 正負をそのまま保持（色分け用）
  
  # ── 5. tidygraph に変換して描画 ──────────────────────────────────────────
  tg <- as_tbl_graph(g)
  
  p <- ggraph(tg, layout = lay) +
    geom_edge_link(
      aes(color = est_sign, width = est_width),
      alpha = edge_alpha
    ) +
    scale_edge_color_gradient2(
      low      = "#3182BD",   # 負 → 青
      mid      = "grey80",
      high     = "#E6550D",   # 正 → 赤
      midpoint = 0,
      name     = "Estimate"
    ) +
    scale_edge_width_identity(guide = "none")  # 幅は凡例不要
  
  if (color_type == "continuous") {
    p <- p +
      geom_node_point(aes(color = color_val, size = node_size)) +
      scale_color_viridis_c(
        option   = ifelse(palette == "viridis", "D", palette),
        na.value = "grey80",
        name     = color_label
      )
  } else {
    p <- p +
      geom_node_point(aes(color = color_cat, size = node_size)) +
      {
        if (is.character(palette) && length(palette) > 1) {
          scale_color_manual(values = palette, name = color_label, na.value = "grey80")
        } else {
          scale_color_brewer(palette = palette, name = color_label, na.value = "grey80")
        }
      }
  }
  
  p <- p +
    scale_size_identity(guide = "legend", name = "Centrality",
                        breaks = range(cent_scaled),
                        labels = round(range(cent), 3))
  
  # ラベル描画（NA のノードは自動的にスキップされる）
  if (!identical(show_labels, "none")) {
    p <- p +
      geom_node_text(
        aes(label = label_name),
        size      = label_size,
        color     = label_color,
        repel     = TRUE,   # ggrepel で重なりを回避
        na.rm     = TRUE
      )
  }
  
  p <- p +
    labs(title = title) +
    theme_void(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
  
  return(p)
}


#!!!!!!!!!!!!!!!!!!!!!
# 使用例
#!!!!!!!!!!!!!!!!!!!!!

## ── 例1: 連続値カラー（log2FC など）、特定遺伝子にラベル ─────────────────
# gene_names <- V(NetworkGraph(grn, graph="umap_graph"))$name
#
# log2fc <- setNames(rnorm(length(gene_names)), gene_names)  # ← 実データに置換
#
# p1 <- plot_grn_custom(
#   grn,
#   node_values  = log2fc,
#   palette      = "RdBu",
#   color_label  = "log2 Fold Change",
#   title        = "GRN UMAP | log2FC",
#   size_range   = c(1, 8),                        # centralityによるサイズ範囲
#   show_labels  = c("GATA1", "TAL1", "MYC")
# )
# print(p1)


## ── 例2: カテゴリカラー、全遺伝子ラベル表示 ─────────────────────────────
# gene_names <- V(NetworkGraph(grn, graph="umap_graph"))$name
# tfs <- GetTFs(grn)
#
# node_cat <- setNames(ifelse(gene_names %in% tfs, "TF", "Target"), gene_names)
#
# p2 <- plot_grn_custom(
#   grn,
#   node_color   = node_cat,
#   palette      = c(TF = "#E64B35", Target = "#4DBBD5"),
#   color_label  = "Gene type",
#   title        = "GRN UMAP | TF vs Target",
#   size_range   = c(1, 10),
#   show_labels  = "all",
#   label_size   = 2.5
# )
# print(p2)


## ── 例3: ラベルなし（デフォルト） ────────────────────────────────────────
# p3 <- plot_grn_custom(
#   grn,
#   node_values  = module_score,
#   palette      = "plasma",
#   color_label  = "Regulon Score",
#   show_labels  = "none"    # デフォルト、ラベルなし
# )
# print(p3)

### EXTRACT Top TF gene 
###
### @input : LINGERres_trans_edges_AddGMMJ , numGene
### @return : TopTFgeme list 
getCorTFgene_pando<- function( pandores  , numGene = 30){
  pandores %>% 
    select(tf, n_genes) %>% 
    distinct() %>% 
    arrange( -n_genes) %>% 
    head(n=numGene) %>%
    pull(TS) %>% 
    return()    

}

getCorTFgene<- function(pandores , numGene = 10 , limit =0 ,coldt = c("red","blue")){
  # pando top
  pandores %>% 
    select(TS, n_genes) %>% 
    distinct() %>% 
    filter(n_genes  >= limit ) %>%
    arrange( -n_genes) %>% 
    head(n=numGene) %>%
    select(TS) %>% 
    #rename( "TS" = "tf")  %>%
    mutate( col = coldt[2], type = "Pando") %>%
    return()    
  
}

### Plot_binHeatmap
###
### @input : heatmapCommon_param,
### gene_order,
### binmatrix , 
### numHighlight=15, 
### numHighlightlimit= 10 , 
### mainname = "",
### # extraparm 
### he_height = unit(4, "cm"),
### ht_width  = unit(8, "cm") ,
### ra_fsize = 10,
### ra_padding = unit(2, "mm"),
### cluster_rows =FALSE 
### @return : SeuratOBJ 

Plot_binHeatmap <- function( heatmapCommon_param,
                             gene_order,
                             binmatrix , 
                             numHighlight=15, 
                             numHighlightlimit= 10 , 
                             mainname = "",
                             # extraparm 
                             he_height = unit(4, "cm"),
                             ht_width  = unit(8, "cm") ,
                             ra_fsize = 10,
                             ra_padding = unit(2, "mm"),
                             cluster_rows =FALSE 
) {
  
  bin_anno_data = heatmapCommon_param$bin_anno_data
  cluster_pallete = heatmapCommon_param$cluster_pallete
  node_info = heatmapCommon_param$node_info
  GeneCluster_col = heatmapCommon_param$GeneCluster_col 
  ManualClusterPalette= heatmapCommon_param$ManualClusterPalette
  HighlightGene = heatmapCommon_param$HighlightGene
  # LINGERres_trans_edges_AddGMM = heatmapCommon_param$LINGERres_trans_edges_AddGMM
  pandores = heatmapCommon_param$pandores  %>% dplyr::rename("TS" ="tf")
  coldt  = heatmapCommon_param$Highlightcoldt
  binmatrix = binmatrix[gene_order,]
  
  ha <- HeatmapAnnotation(
    `  PredictedDay` = anno_barplot(
      bin_anno_data,
      bar_width = 1,
      gp = gpar(fill = cluster_pallete[colnames(bin_anno_data)] ),
      border = TRUE,
      axis = TRUE,
      height = he_height
    )
  )
  row_ha = rowAnnotation( GeneCluster = node_info[gene_order,"Cluster"] ,
                          #SubCluster = node_info[gene_order,"SubCluster"] ,
                          #ManualCluster = node_info[gene_order,"ManualCluster"] ,
                          col = list(GeneCluster = GeneCluster_col 
                                    # SubCluster = GeneCluster_col ,
                                    # ManualCluster = ManualClusterPalette 
                          ) )
  

  PnadoTF <- getCorTFgene(pandores %>% filter(TS %in% gene_order )  , 
               numGene=numHighlight, limit= numHighlightlimit , coldt = coldt)
  #print(PnadoTF)
  plotHighlightGene <- c( HighlightGene , PnadoTF)
  
  #print(plotHighlightGene)
  HighlightGene_col <- plotHighlightGene$col
  names(HighlightGene_col) <- plotHighlightGene$TS
  HighlightGene_locus <- which(gene_order %in% plotHighlightGene$TS)
  
  # 転写因子の上位を作る必要あり
  ht=Heatmap(binmatrix %>% as.matrix() ,
             name = mainname, 
             top_annotation = ha, 
             left_annotation = row_ha,
             show_column_names = F,
             show_row_names =  F,
             cluster_columns = FALSE,
             cluster_rows = cluster_rows , width = ht_width )+
    rowAnnotation(link = anno_mark(at = HighlightGene_locus  , 
                                   labels = gene_order[HighlightGene_locus ], 
                                   labels_gp = gpar(col = HighlightGene_col[gene_order[HighlightGene_locus ]] 
                                                    ,fontsize = ra_fsize), padding = ra_padding ) )

  return(ht)
}



