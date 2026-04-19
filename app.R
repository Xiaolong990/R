# ======================
# Shiny App: 膜过滤数据数值微分分析平台 (科普增强版)
# ======================
library(shiny)
library(readxl)
library(dplyr)
library(plotly)
library(splines)
library(signal)
library(pracma)
library(numDeriv)
library(reshape2)
library(bslib) # 用于更好的主题和工具提示

# ----------------------
# 1. 知识库：数学原理与参数说明
# ----------------------
method_info <- list(
  center_diff = list(
    title = "中心差分法 (Centered Difference)",
    math = "$$ f'(x_i) \\approx \\frac{f(x_{i+h}) - f(x_{i-h})}{x_{i+h} - x_{i-h}} $$",
    desc = "最基础的数值微分方法。利用当前点两侧距离为 $h$ 的点来估算斜率。
    <br><strong>优点</strong>：计算极快，逻辑简单，对线性区域非常准确。
    <br><strong>缺点</strong>：对噪声极度敏感（噪声会被放大），分辨率受 $h$ 限制。",
    usage = "适用于数据非常干净（低噪声）的情况，或者作为其他高级方法的基准对比。"
  ),
  spline = list(
    title = "平滑样条 (Smoothing Spline)",
    math = "最小化目标函数：$$\\sum (y_i - f(x_i))^2 + \\lambda \\int (f''(t))^2 dt$$",
    desc = "通过构建一条光滑的分段多项式曲线来拟合数据，然后对这条曲线进行解析求导。
    <br><strong>优点</strong>：生成的导数曲线非常光滑，全局连续性好。
    <br><strong>缺点</strong>：计算量较大，边界处可能出现伪影（Runge 现象的变体）。",
    usage = "适用于需要极其光滑的导数曲线，且数据存在中等程度噪声的情况。"
  ),
  loess = list(
    title = "Loess 局部回归 (Locally Estimated Scatterplot Smoothing)",
    math = "在每个点 $$x_0$$ 附近，利用加权最小二乘法拟合一个低阶多项式。",
    desc = "一种非参数回归方法。它不假设全局函数形式，而是在每个点的邻域内局部拟合。
    <br><strong>优点</strong>：能很好地适应数据的局部特征和非线性变化，鲁棒性强。
    <br><strong>缺点</strong>：计算速度慢（大数据集慎用），边界效应较明显。",
    usage = "适用于数据具有复杂的局部波动，且噪声分布不均匀的情况。"
  ),
  sg = list(
    title = "Savitzky-Golay 滤波器",
    math = "在滑动窗口内拟合多项式：$$y = a_0 + a_1 x + \\dots + a_p x^p$$，导数即 $$a_1$$。",
    desc = "本质是在滑动窗口内进行多项式最小二乘拟合。它能保留信号的高频特征（如峰宽、峰高）。
    <br><strong>优点</strong>：在平滑噪声的同时，能很好地保留信号的形状特征（矩）。
    <br><strong>缺点</strong>：严格要求数据点等间距（或近似等间距），否则误差巨大。",
    usage = "适用于色谱、光谱或膜过滤等近似等间距采样且需要保留峰值特征的数据。"
  )
)

param_info <- list(
  h = "差分步长 (h)。决定中心差分取点的跨度。<br>• <strong>调大</strong>：抗噪能力增强，但会模糊细节，降低时间分辨率。<br>• <strong>调小</strong>：分辨率高，但对噪声极其敏感，曲线会剧烈震荡。",
  sigma = "高斯平滑系数 (σ)。对差分后的结果进行后处理平滑。<br>• <strong>调大</strong>：曲线更平滑，但可能滞后或抹平尖峰。<br>• <strong>调小</strong>：保留更多细节，但残留噪声多。",
  spar = "平滑参数 (spar, 0~1)。控制样条曲线的拟合程度。<br>• <strong>调大 (接近1)</strong>：曲线更直、更平滑，可能欠拟合。<br>• <strong>调小 (接近0)</strong>：曲线穿过更多数据点，可能过拟合噪声。",
  span = "局部回归跨度 (Span, 0~1)。决定 Loess 拟合时参与计算的数据比例。<br>• <strong>调大</strong>：考虑更多点，曲线更平滑，反应迟钝。<br>• <strong>调小</strong>：只考虑极少数邻近点，能捕捉快速变化，但易受噪声干扰。",
  sg_p = "多项式阶数 (p)。SG 滤波器拟合的多项式次数。<br>• 通常选 2 (二次) 或 4 (四次)。阶数过高会导致过拟合噪声。",
  sg_n = "窗口大小 (n)。SG 滤波器滑动窗口的点数。<br>• <strong>调大</strong>：平滑效果强，但可能失真。<br>• <strong>注意</strong>：必须是奇数，且应小于数据总量的 10%-20%。"
)

# ----------------------
# 2. 辅助函数
# ----------------------

my_gauss_filt <- function(x, sigma = 2.0) {
  if (sigma <= 0) return(x)
  kernel_size <- ceiling(6 * sigma)
  if (kernel_size %% 2 == 0) kernel_size <- kernel_size + 1
  if (kernel_size < 3) kernel_size <- 3
  half_win <- (kernel_size - 1) / 2
  x_kernel <- seq(-half_win, half_win)
  weights <- dnorm(x_kernel, mean = 0, sd = sigma)
  weights <- weights / sum(weights)
  filtered_x <- stats::filter(x, filter = weights, method = "convolution", sides = 2)
  return(filtered_x)
}

process_data <- function(file_path) {
  if (is.null(file_path)) return(NULL)
  tryCatch({
    raw_data <- read_excel(file_path$datapath, col_names = TRUE)
    if (nrow(raw_data) == 0) stop("Excel 文件为空。")
    original_names <- names(raw_data)
    clean_names <- tolower(trimws(original_names))
    t_keywords <- c("t", "time", "min", "hour", "sec", "second")
    v_keywords <- c("v", "vol", "volume", "liter", "l", "m3")
    t_idx <- NULL; v_idx <- NULL
    for (i in seq_along(clean_names)) {
      name <- clean_names[i]
      if (is.null(t_idx) && any(sapply(t_keywords, function(k) grepl(k, name)))) t_idx <- i
      if (is.null(v_idx) && any(sapply(v_keywords, function(k) grepl(k, name)))) v_idx <- i
    }
    if (is.null(t_idx) || is.null(v_idx)) {
      is_num_col <- sapply(raw_data, function(x) sum(!is.na(suppressWarnings(as.numeric(x)))) > length(x) * 0.8)
      numeric_indices <- which(is_num_col)
      if (length(numeric_indices) >= 2) {
        if (is.null(t_idx)) t_idx <- numeric_indices[1]
        if (is.null(v_idx)) v_idx <- numeric_indices[2]
      } else stop("无法自动识别列。")
    }
    col_t <- original_names[t_idx]; col_v <- original_names[v_idx]
    clean_data <- data.frame(
      t = suppressWarnings(as.numeric(raw_data[[col_t]])),
      V = suppressWarnings(as.numeric(raw_data[[col_v]]))
    ) 
    if (nrow(clean_data) < 5) stop("数据点太少。")
    clean_data$t_over_V <- clean_data$t / clean_data$V
    return(list(data = clean_data, col_t_name = col_t, col_v_name = col_v))
  }, error = function(e) { stop(paste("读取失败:", e$message)) })
}

compute_all_derivatives <- function(data, params) {
  V <- data$V; t <- data$t; y <- data$t_over_V; n <- length(V)
  res <- data.frame(t = t, V = V, t_over_V = y)
  h <- params$diff_h
  if (h < 1) h <- 1
  if (h >= n/2) h <- floor((n-1)/2)
  
  # Center Diff
  cd_raw <- rep(NA, n)
  if (n > 2*h) {
    for (i in (h+1):(n-h)) {
      denom <- V[i+h] - V[i-h]
      if (denom != 0) cd_raw[i] <- (y[i+h] - y[i-h]) / denom
    }
  }
  for (j in 1:h) { if (j+1 <= n && (V[j+1]-V[j])!=0) cd_raw[j] <- (y[j+1] - y[j]) / (V[j+1] - V[j]) }
  for (k in (n-h+1):n) { if (k-1 >= 1 && (V[k]-V[k-1])!=0) cd_raw[k] <- (y[k] - y[k-1]) / (V[k] - V[k-1]) }
  res$Center_Diff <- my_gauss_filt(cd_raw, sigma = params$gauss_sigma)
  
  # Spline
  tryCatch({
    spline_fit <- smooth.spline(V, y, spar = params$spline_spar)
    res$Smooth_Spline <- predict(spline_fit, V, deriv = 1)$y
  }, error = function(e) { res$Smooth_Spline <- NA })
  
  # Loess
  tryCatch({
    safe_span <- max(min(params$loess_span, 0.99), 2/n + 0.01)
    loess_model <- loess(y ~ V, span = safe_span, degree = 2)
    loess_func <- function(x) predict(loess_model, newdata = data.frame(V = x))
    res$Loess <- sapply(V, function(x) tryCatch(grad(loess_func, x, method = "Richardson"), error = function(e) NA))
  }, error = function(e) { res$Loess <- NA })
  
  # SG
  res$SG <- NA
  if (n > 5) {
    dV <- diff(V)
    if (min(dV) > 0 && (max(dV)/min(dV)) < 1.5) {
      tryCatch({
        win_size <- params$sg_window
        if (win_size %% 2 == 0) win_size <- win_size + 1
        win_size <- min(win_size, n)
        if (win_size <= params$sg_poly) win_size <- params$sg_poly + 1 + (params$sg_poly + 1) %% 2
        res$SG <- sgolayfilt(y, p = params$sg_poly, n = win_size, m = 1, ts = mean(dV))
      }, error = function(e) { res$SG <- NA })
    }
  }
  return(res)
}

generate_plots <- function(data, deriv_df) {
  plots <- list()
  plots$p_original <- plot_ly(data, x = ~t, y = ~t_over_V, type = 'scatter', mode = 'markers+lines', name = 't/V', marker = list(size = 6), line = list(color = '#0072B2')) %>%
    layout(title = '原始数据 (t/V vs t)', xaxis = list(title = 'Time (t)'), yaxis = list(title = 't/V'), margin = list(l=50, r=50, t=50, b=50))
  
  deriv_long <- melt(deriv_df, id.vars = c("t", "V"), variable.name = "Method", value.name = "Value", na.rm = TRUE)
  cols <- c("Center_Diff" = "#D62728", "Smooth_Spline" = "#1F77B4", "Loess" = "#2CA02C", "SG" = "#9467BD")
  
  p_compare <- plot_ly()
  for (m in names(cols)) {
    sub <- deriv_long[deriv_long$Method == m, ]
    if (nrow(sub) > 0 && !all(is.na(sub$Value))) {
      p_compare <- p_compare %>% add_trace(data = sub, x = ~t, y = ~Value, type = 'scatter', mode = 'lines', name = m, line = list(color = cols[m], width = 2), inherit = FALSE)
    }
  }
  plots$p_compare <- p_compare %>% layout(title = '四种微分方法结果对比 (vs Time)', xaxis = list(title = 'Time (t)'), yaxis = list(title = 'd(t/V)/dt'), hovermode = 'closest', legend = list(orientation = 'h'), margin = list(l=50, r=50, t=50, b=50))
  
  plots$p_cd <- plot_ly(deriv_df, x = ~t, y = ~Center_Diff, type = 'scatter', mode = 'lines', name = 'Center Diff', line = list(color = '#D62728')) %>%
    layout(title = '中心差分详情', xaxis = list(title = 'Time (t)'), yaxis = list(title = 'Derivative'), margin = list(l=50, r=50, t=50, b=50))
  
  other_methods <- deriv_long[deriv_long$Method != "Center_Diff", ]
  plots$p_others <- plot_ly()
  for (m in unique(other_methods$Method)) {
    sub <- other_methods[other_methods$Method == m, ]
    if (nrow(sub) > 0 && !all(is.na(sub$Value))) {
      plots$p_others <- plots$p_others %>% add_trace(data = sub, x = ~t, y = ~Value, type = 'scatter', mode = 'lines', name = m, line = list(color = cols[m]), inherit = FALSE)
    }
  }
  plots$p_others <- plots$p_others %>% layout(title = '高级算法对比', xaxis = list(title = 'Time (t)'), yaxis = list(title = 'Derivative'), margin = list(l=50, r=50, t=50, b=50))
  return(plots)
}

# ----------------------
# 3. UI 界面设计
# ----------------------
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  titlePanel("膜过滤数据数值微分分析平台 (科普版)"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("1. 数据导入"),
      fileInput("file_excel", "上传 Excel (.xlsx)", accept = ".xlsx"),
      hr(),
      
      h4("2. 参数设置"),
      p("点击 '?' 查看数学原理与参数影响说明。"),
      
      wellPanel(
        h5(
          "中心差分 (Centered Diff)", 
          actionButton("help_center", "?", style = "padding: 2px 8px; font-size: 12px; float: right;")
        ),
        sliderInput("param_h", tags$span("差分步长 ", bslib::tooltip("h", param_info$h)), min = 1, max = 20, value = 2, step = 1),
        sliderInput("param_gauss", tags$span("高斯平滑 σ", bslib::tooltip("sigma", param_info$sigma)), min = 0, max = 10, value = 1.0, step = 0.5)
      ),
      
      wellPanel(
        h5(
          "平滑样条 (Spline)", 
          actionButton("help_spline", "?", style = "padding: 2px 8px; font-size: 12px; float: right;")
        ),
        sliderInput("param_spline", tags$span("平滑参数 spar", bslib::tooltip("spar", param_info$spar)), min = 0, max = 1, value = 0.5, step = 0.05)
      ),
      
      wellPanel(
        h5(
          "Loess 局部回归", 
          actionButton("help_loess", "?", style = "padding: 2px 8px; font-size: 12px; float: right;")
        ),
        sliderInput("param_loess", tags$span("跨度 Span", bslib::tooltip("span", param_info$span)), min = 0.1, max = 0.9, value = 0.3, step = 0.05)
      ),
      
      wellPanel(
        h5(
          "Savitzky-Golay", 
          actionButton("help_sg", "?", style = "padding: 2px 8px; font-size: 12px; float: right;")
        ),
        sliderInput("param_sg_p", tags$span("多项式阶数 p", bslib::tooltip("p", param_info$sg_p)), min = 1, max = 4, value = 2, step = 1),
        sliderInput("param_sg_n", tags$span("窗口大小 n", bslib::tooltip("n", param_info$sg_n)), min = 3, max = 21, value = 7, step = 2)
      ),
      
      hr(),
      actionButton("btn_calc", "重新计算", class = "btn-primary", style = "width: 100%"),
      br(), br(),
      h4("3. 数据导出"),
      downloadButton("download_data", "下载计算结果 (CSV)", class = "btn-success", style = "width: 100%"),
      br(), br(),
      verbatimTextOutput("info_status")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("📈 综合对比视图", 
                 fluidRow(column(12, div(style = "border: 1px solid #ddd; padding: 10px; border-radius: 5px; background-color: #fff;", plotlyOutput("plot_compare", height = "550px")))),
                 fluidRow(
                   column(6, div(style = "border: 1px solid #ddd; padding: 10px; border-radius: 5px; background-color: #fff;", plotlyOutput("plot_original", height = "400px"))),
                   column(6, div(style = "border: 1px solid #ddd; padding: 10px; border-radius: 5px; background-color: #fff;", plotlyOutput("plot_others", height = "400px")))
                 )
        ),
        tabPanel("🔍 方法细节视图",
                 fluidRow(
                   column(6, div(style = "border: 1px solid #ddd; padding: 10px; border-radius: 5px; background-color: #fff;", plotlyOutput("plot_cd", height = "500px"))),
                   column(6, div(style = "border: 1px solid #ddd; padding: 10px; border-radius: 5px; background-color: #fff;", plotlyOutput("plot_others_detail", height = "500px")))
                 )
        ),
        # 新增：原理科普标签页
        tabPanel("📚 方法原理与指南",
                 fluidRow(
                   column(12,
                          h3("数值微分方法原理解析"),
                          p("本工具提供了四种主流的数值微分算法。以下是它们的数学原理、优缺点及适用场景。"),
                          hr()
                   )
                 ),
                 fluidRow(
                   column(6,
                          div(style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; border-left: 5px solid #D62728;",
                              h4(method_info$center_diff$title, style = "color: #D62728;"),
                              withMathJax(helpText(method_info$center_diff$math)),
                              p(HTML(method_info$center_diff$desc)),
                              p(HTML(method_info$center_diff$usage))
                          )
                   ),
                   column(6,
                          div(style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; border-left: 5px solid #1F77B4;",
                              h4(method_info$spline$title, style = "color: #1F77B4;"),
                              withMathJax(helpText(method_info$spline$math)),
                              p(HTML(method_info$spline$desc)),
                              p(HTML(method_info$spline$usage))
                          )
                   )
                 ),
                 fluidRow(
                   column(6,
                          div(style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; border-left: 5px solid #2CA02C;",
                              h4(method_info$loess$title, style = "color: #2CA02C;"),
                              withMathJax(helpText(method_info$loess$math)),
                              p(HTML(method_info$loess$desc)),
                              p(HTML(method_info$loess$usage))
                          )
                   ),
                   column(6,
                          div(style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; border-left: 5px solid #9467BD;",
                              h4(method_info$sg$title, style = "color: #9467BD;"),
                              withMathJax(helpText(method_info$sg$math)),
                              p(HTML(method_info$sg$desc)),
                              p(HTML(method_info$sg$usage))
                          )
                   )
                 )
        )
      )
    )
  )
)

# ----------------------
# 4. 服务器逻辑
# ----------------------
server <- function(input, output, session) {
  data_clean <- reactiveVal(NULL)
  data_meta <- reactiveVal(NULL)
  deriv_results <- reactiveVal(NULL)
  
  # 模态框逻辑：显示详细说明
  observeEvent(input$help_center, {
    showModal(modalDialog(
      title = method_info$center_diff$title,
      withMathJax(helpText(method_info$center_diff$math)),
      p(HTML(method_info$center_diff$desc)),
      p(HTML(param_info$h)),
      p(HTML(param_info$sigma)),
      easyClose = TRUE, footer = NULL
    ))
  })
  
  observeEvent(input$help_spline, {
    showModal(modalDialog(
      title = method_info$spline$title,
      withMathJax(helpText(method_info$spline$math)),
      p(HTML(method_info$spline$desc)),
      p(HTML(param_info$spar)),
      easyClose = TRUE, footer = NULL
    ))
  })
  
  observeEvent(input$help_loess, {
    showModal(modalDialog(
      title = method_info$loess$title,
      withMathJax(helpText(method_info$loess$math)),
      p(HTML(method_info$loess$desc)),
      p(HTML(param_info$span)),
      easyClose = TRUE, footer = NULL
    ))
  })
  
  observeEvent(input$help_sg, {
    showModal(modalDialog(
      title = method_info$sg$title,
      withMathJax(helpText(method_info$sg$math)),
      p(HTML(method_info$sg$desc)),
      p(HTML(param_info$sg_p)),
      p(HTML(param_info$sg_n)),
      easyClose = TRUE, footer = NULL
    ))
  })
  
  observeEvent(input$file_excel, {
    req(input$file_excel)
    showNotification("正在读取文件...", type = "message", duration = 2)
    tryCatch({
      result_list <- process_data(input$file_excel)
      data_clean(result_list$data)
      data_meta(list(t_name = result_list$col_t_name, v_name = result_list$col_v_name))
      showNotification(paste("成功读取:", result_list$col_t_name, "&", result_list$col_v_name), type = "message") 
    }, error = function(e) {
      showNotification(paste("错误:", e$message), type = "error", duration = 10)
      data_clean(NULL); data_meta(NULL); deriv_results(NULL)
    })
  })
  
  observeEvent({input$btn_calc; data_clean()}, {
    req(data_clean())
    isolate({
      params <- list(
        diff_h = input$param_h,       
        gauss_sigma = input$param_gauss,
        spline_spar = input$param_spline,
        loess_span = input$param_loess,
        sg_poly = input$param_sg_p,
        sg_window = input$param_sg_n
      )
      tryCatch({
        res <- compute_all_derivatives(data_clean(), params)
        deriv_results(res)
        showNotification("计算完成", type = "message", duration = 1)
      }, error = function(e) {
        showNotification(paste("计算错误:", e$message), type = "error")
        deriv_results(NULL)
      })
    })
  }, ignoreInit = FALSE)
  
  output$info_status <- renderText({
    if (is.null(data_clean())) return("等待上传文件...")
    meta <- data_meta()
    n_pts <- nrow(data_clean())
    if (!is.null(meta)) sprintf("列: %s, %s\n数据点: %d\n状态: 就绪", meta$t_name, meta$v_name, n_pts)
    else sprintf("数据点: %d\n状态: 就绪", n_pts)
  })
  
  output$plot_original <- renderPlotly({ req(deriv_results()); req(data_clean()); generate_plots(data_clean(), deriv_results())$p_original })
  output$plot_compare <- renderPlotly({ req(deriv_results()); req(data_clean()); generate_plots(data_clean(), deriv_results())$p_compare })
  output$plot_cd <- renderPlotly({ req(deriv_results()); req(data_clean()); generate_plots(data_clean(), deriv_results())$p_cd })
  output$plot_others <- renderPlotly({ req(deriv_results()); req(data_clean()); generate_plots(data_clean(), deriv_results())$p_others })
  output$plot_others_detail <- renderPlotly({ req(deriv_results()); req(data_clean()); generate_plots(data_clean(), deriv_results())$p_others })
  
  output$download_data <- downloadHandler(
    filename = function() { paste0("membrane_diff_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv") },
    content = function(file) { req(deriv_results()); write.csv(deriv_results(), file, row.names = FALSE) }
  )
}

shinyApp(ui = ui, server = server)