use crate::global_variables::*;

pub struct PostResult {
    pub name: String,
    pub label: String,
    pub value: Float,
    pub unit: Option<String>,
}

impl PostResult {
    pub fn new(name: String, label: String, value: Float, unit: Option<String>) -> Self {
        Self {
            name,
            label,
            value,
            unit,
        }
    }
}
