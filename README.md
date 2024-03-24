
## Features

☢ Deeply Reactive, Directly Mutate State at any level to Update Component

🌿 Always Fresh State

☕ Zero Dependencies, Ultra Light-Weight `< 1kb`

<br/>




## 🌻 Motivation

This is a collection of script in Python, the example is **DCDC state space matrix for buck converter**.
This project allows you to have immediate physical feedback on reality when the system parameters change 


## Installation
No installation needed.

```bash

```
<br/>



## 🌿 State is always fresh !

Start always a fresh windows, this can help to avoid some issues.

#### What does that mean ?




## ❓ FAQs

<!-- faq 1 -->
<details>
<summary>Can I use useRS hook more than once ? </summary>
<br/>

**Yes.** You don't have to put all of the state of the component inside the state object. You can use the hook more than once.

#### Example

```javascript
const todos = useRS([])
const form = useRS({
  name: '',
  age: 0,
})
```

While this is okay, **I would advise you to not do this**, Because putting all of state in one object gives you **better *performance** in the case of radioactive-state. (because of better mutation batching)

It would also be **hard to store simple value types**, because simple value types can not be mutated and so you would need to wrap it inside an object.

#### Example

```javascript
const count = useRS(0) // invalid, gives error ❌
```

```javascript
const count = useRS( { value: 0 }) // works ✅
```

This would also make creating reactive bindings awkward. That's why it is **strongly recommended to store all the state into a single object** by using useRS only once !

---
</details>


<!-- FAQ 2 -->
<details>
<summary> Is this magic, How does it work ? </summary>
<br/>
radioactive-state uses **JavaScript Proxy** to create a deeply reactive state by recursively proxifying the state. Whenever a mutation occurs in the state tree, a function is called with information about where the mutation took place which schedules an async re-render to update the component to reflect the changes in state to UI.
</details>
<br/>




## 💙 Contributing

PR's are welcome !

Found a Bug ? Create an Issue.

Chat on [Linkedin](https://www.linkedin.com/in/roberto-di-lorenzo-phd-0b841997/ "Roberto Di Lorenzo")

<br/>




## 💖 Like this project ?

Leave a ⭐ If you think this project is cool.

[Share with the world](https://twitter.com/intent/tweet?url=https%3A%2F%2Fgithub.com%2FMananTank%2Fradioactive-state&via=MananTank_&text=Make%20your%20@react%20App%20Truly%20Reactive%20with%20radioactive-state&hashtags=react%2CradioactiveState) ✨

<br/>




## 👨‍💻 Author

### Roberto Di Lorenzo

[Linkedin](https://www.linkedin.com/in/roberto-di-lorenzo-phd-0b841997/ "Roberto Di Lorenzo")

<br/>




## 🍁 Licence

**ISC**
