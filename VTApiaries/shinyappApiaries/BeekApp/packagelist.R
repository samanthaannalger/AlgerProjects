package_df<-as.data.frame(installed.packages("/Library/Frameworks/R.framework/Versions/3.3/Resources/library"))

package_list <- as.character(package_df$Package)
write.csv(package_list, file = "packagelist.csv")

install.packages(package_list)