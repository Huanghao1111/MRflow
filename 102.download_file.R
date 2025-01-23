#You can use any download tool(IDM aria2 .etc) to download the complete data,
# as it is much faster than using the code below.
download_data <- function(file_name, output_path){
  if(!dir.exists(output_path)){
    dir.create(output_path, recursive = T)
  }
  urls <- read.table(file_name)
  for(url in urls$V1){
    a <- basename(url)
    message("download file:", a)
    curl::curl_download(url, file.path(output_path, a), quiet = FALSE)
  }
}
# exposure
download_data("url_txt/1400met_urls.txt", "data/01origin_data/1400met")
download_data("url_txt/233met_urls.txt", "data/01origin_data/233met")
download_data("url_txt/338met_urls.txt", "data/01origin_data/338met")

#outcome

