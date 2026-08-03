# 数据字典

公开仓库中的 CSV 是整理后的美国月度宏观数据。最终模型只使用下列 5 个序列；文件中的其他列为早期研究阶段收集的辅助变量，当前展示版本没有使用。

| CSV 字段 | FRED 序列 | 含义 | 原始单位 | 模型中的处理方式 |
|---|---|---|---|---|
| `Income` | [DSPI](https://fred.stlouisfed.org/series/DSPI) | 个人可支配收入 | 十亿美元，季调年率 | `log(DSPI / PCEPI x 100)` |
| `Price` | [PCEPI](https://fred.stlouisfed.org/series/PCEPI) | 个人消费支出价格指数 | 指数，2017 年 = 100 | 作为价格平减指数 |
| `Consumption` | [PCE](https://fred.stlouisfed.org/series/PCE) | 个人消费支出 | 十亿美元，季调年率 | `log(PCE / PCEPI x 100)` |
| `Real_income_excptTransfer` | [W875RX1](https://fred.stlouisfed.org/series/W875RX1) | 扣除经常转移收入后的实际个人收入 | 2017 年不变价十亿美元，季调年率 | 取自然对数 |
| `Unemploy_rate` | [UNRATE](https://fred.stlouisfed.org/series/UNRATE) | 美国民用失业率 | 百分比，季调 | 除以 100 后提取失业率缺口 |

## 转移支付因子的构造

模型使用以下对数差值近似衡量税收与转移支付对居民可支配资源的相对贡献：

```text
transfer_factor = log(实际可支配收入)
                - log(扣除转移支付后的实际个人收入)
```

转移支付会随失业和经济周期自动变化，因此模型中的对应系数只表示控制其他变量后的相关关系，不能直接理解为因果财政乘数。

## 样本与数据说明

- 原始文件覆盖区间：1959 年 1 月至 2025 年 9 月；
- 归档模型使用区间：前 733 个月，即 1959 年 1 月至 2020 年 1 月；
- 数据频率：月度；
- 原始研究流程没有记录准确的下载日期和 FRED 历史版本；
- 部分未使用的辅助变量存在缺失值，不影响当前模型使用的 5 个核心序列。
