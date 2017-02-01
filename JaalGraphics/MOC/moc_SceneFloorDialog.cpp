/****************************************************************************
** Meta object code from reading C++ file 'SceneFloorDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/SceneFloorDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'SceneFloorDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JSceneFloorDialog_t {
    QByteArrayData data[12];
    char stringdata0[123];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JSceneFloorDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JSceneFloorDialog_t qt_meta_stringdata_JSceneFloorDialog = {
    {
QT_MOC_LITERAL(0, 0, 17), // "JSceneFloorDialog"
QT_MOC_LITERAL(1, 18, 10), // "setPattern"
QT_MOC_LITERAL(2, 29, 0), // ""
QT_MOC_LITERAL(3, 30, 8), // "setColor"
QT_MOC_LITERAL(4, 39, 12), // "checkDisplay"
QT_MOC_LITERAL(5, 52, 12), // "setDirection"
QT_MOC_LITERAL(6, 65, 11), // "setDistance"
QT_MOC_LITERAL(7, 77, 9), // "setLength"
QT_MOC_LITERAL(8, 87, 8), // "setLines"
QT_MOC_LITERAL(9, 96, 13), // "keyPressEvent"
QT_MOC_LITERAL(10, 110, 10), // "QKeyEvent*"
QT_MOC_LITERAL(11, 121, 1) // "e"

    },
    "JSceneFloorDialog\0setPattern\0\0setColor\0"
    "checkDisplay\0setDirection\0setDistance\0"
    "setLength\0setLines\0keyPressEvent\0"
    "QKeyEvent*\0e"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JSceneFloorDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   54,    2, 0x08 /* Private */,
       3,    0,   55,    2, 0x08 /* Private */,
       4,    0,   56,    2, 0x08 /* Private */,
       5,    0,   57,    2, 0x08 /* Private */,
       6,    0,   58,    2, 0x08 /* Private */,
       7,    0,   59,    2, 0x08 /* Private */,
       8,    0,   60,    2, 0x08 /* Private */,
       9,    1,   61,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 10,   11,

       0        // eod
};

void JSceneFloorDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JSceneFloorDialog *_t = static_cast<JSceneFloorDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->setPattern(); break;
        case 1: _t->setColor(); break;
        case 2: _t->checkDisplay(); break;
        case 3: _t->setDirection(); break;
        case 4: _t->setDistance(); break;
        case 5: _t->setLength(); break;
        case 6: _t->setLines(); break;
        case 7: _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObject JSceneFloorDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JSceneFloorDialog.data,
      qt_meta_data_JSceneFloorDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JSceneFloorDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JSceneFloorDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JSceneFloorDialog.stringdata0))
        return static_cast<void*>(const_cast< JSceneFloorDialog*>(this));
    if (!strcmp(_clname, "Ui::SceneFloorDialog"))
        return static_cast< Ui::SceneFloorDialog*>(const_cast< JSceneFloorDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JSceneFloorDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 8)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 8;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
